#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include <xrt/xrt_bo.h>
#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>

#include "spot_finder.h"

int dispersion_filter(std::vector<uint16_t> &image_in,
                      std::vector<uint16_t> &mask_out) {
  int32_t knl = 3;

  float sigma_b = 6.0f, sigma_s = 3.0f;

  std::vector<uint32_t> m_sat(NY * NX);
  std::vector<uint32_t> i_sat(NY * NX);
  std::vector<uint32_t> i2_sat(NY * NX);

  // copy image data in, fill SATs
  for (int32_t i = 0, k = 0; i < NY; i++) {
    uint32_t _m = 0, _i = 0, _i2 = 0;
    for (int32_t j = 0; j < NX; j++, k++) {
      uint32_t m = image_in[k] > 0xfffd ? 0 : 1;
      uint32_t p = m * image_in[k];
      _m += m;
      _i += p;
      _i2 += p * p;
      m_sat[k] = i > 0 ? _m + m_sat[k - NX] : _m;
      i_sat[k] = i > 0 ? _i + i_sat[k - NX] : _i;
      i2_sat[k] = i > 0 ? _i2 + i2_sat[k - NX] : _i2;
    }
  }

  // roll over now computing the mean, variance etc.
  for (int32_t i = 0, k = 0; i < NY; i++) {
    for (int32_t j = 0; j < NX; j++, k++) {
      int32_t j0 = j - knl - 1;
      int32_t j1 = j < (NX - knl) ? j + knl : NX - 1;
      int32_t i0 = i - knl - 1;
      int32_t i1 = i < (NY - knl) ? i + knl : NY - 1;

      int32_t a = i1 * NX + j1;
      int32_t b = i0 * NX + j1;
      int32_t c = i1 * NX + j0;
      int32_t d = i0 * NX + j0;

      uint32_t m_sum = m_sat[a], i_sum = i_sat[a], i2_sum = i2_sat[a];

      if (j0 >= 0 && i0 >= 0) {
        m_sum += m_sat[d] - m_sat[b] - m_sat[c];
        i_sum += i_sat[d] - i_sat[b] - i_sat[c];
        i2_sum += i2_sat[d] - i2_sat[b] - i2_sat[c];
      } else if (j0 >= 0) {
        m_sum -= m_sat[c];
        i_sum -= i_sat[c];
        i2_sum -= i2_sat[c];
      } else if (i0 >= 0) {
        m_sum -= m_sat[b];
        i_sum -= i_sat[b];
        i2_sum -= i2_sat[b];
      }

      uint16_t signal = 0, isignal = 0;
      uint32_t p = image_in[k] > 0xfffd ? 0 : image_in[k];

      if (p > 0 && m_sum >= 2) {
        double bg_lhs = (double)m_sum * i2_sum - (double)i_sum * i_sum -
                        (double)i_sum * (m_sum - 1);
        double bg_rhs = i_sum * sigma_b * sqrt((double)2.0f * (m_sum - 1));
        uint16_t background = bg_lhs > bg_rhs;
        double fg_lhs = (double)m_sum * p - (double)i_sum;
        double fg_rhs = sigma_s * sqrtf((double)i_sum * m_sum);
        uint16_t foreground = fg_lhs > fg_rhs;
        signal = background && foreground;
      }


      if (p > 0 && m_sum >= 2) {
        int64_t bg_lhs = (int64_t)m_sum * i2_sum - (int64_t)i_sum * i_sum -
                         (int64_t)i_sum * (m_sum - 1);
        int64_t bg_rhs2 = i_sum * i_sum * 72 * (m_sum - 1);
        uint8_t background = (bg_lhs > 0) && ((bg_lhs * bg_lhs) > bg_rhs2);
        int32_t fg_lhs = (int32_t)m_sum * p - (int32_t)i_sum;
        int32_t fg_rhs2 = 9 * i_sum * m_sum;
        uint8_t foreground = (fg_lhs > 0) && ((fg_lhs * fg_lhs) > fg_rhs2);
        isignal = background && foreground;
      }

      if (isignal != signal) {
        std::cout << "Algorithm fail at " << k << std::endl;
      }

      // save the pixel value for later use in connected component labelling
      mask_out[k] = signal * p;
    }
  }

  return 0;
}

void random_image(std::vector<uint16_t> &image, float rate, uint32_t seed) {
  std::mt19937 generator(seed);
  std::poisson_distribution<int> background(rate);
  std::poisson_distribution<int> spot_counts(100 * rate);
  std::normal_distribution<float> spot_xy(4, 1);
  std::uniform_int_distribution<int> bad_k(0, (NX * NY) - 1);

  int NY2 = 32 * (NY / 32);
  int NX2 = 32 * (NX / 32);

  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      image[i * NX + j] = (uint16_t)background(generator);
    }
  }

  for (int i = 0; i < NY2; i += 32) {
    for (int j = 0; j < NX2; j += 32) {
      int counts = spot_counts(generator);
      while (counts) {
        int x = (int)spot_xy(generator);
        int y = (int)spot_xy(generator);
        if (x >= 0 && x < 8 && y >= 0 && y < 8) {
          counts -= 1;
          ++image[(i + y + 12) * NX + j + x + 12];
        }
      }
    }
  }

  for (int b = 0; b < 128; b++) {
    image[bad_k(generator)] = 0xffff;
  }
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << argv[0] << " kernel.xclbin" << std::endl;
    return 1;
  }

  // only consider single device right now
  auto device = xrt::device(0);
  auto uuid = device.load_xclbin(argv[1]);

  auto krnl = xrt::kernel(device, uuid, "spot_finder");

  auto bo_in = xrt::bo(device, NY * NX * sizeof(uint16_t), krnl.group_id(0));
  auto bo_out = xrt::bo(device, NY * NX * sizeof(uint16_t), krnl.group_id(1));

  auto image_in = (uint16_t *)bo_in.map<pack16_t *>();
  auto mask = (uint16_t *)bo_out.map<pack16_t *>();

  std::vector<uint16_t> image(NY * NX);
  std::vector<uint16_t> ref(NY * NX);

  int wrong = 0;

  auto t0 = std::chrono::high_resolution_clock::now();
  random_image(image, 10, 50);
  auto t1 = std::chrono::high_resolution_clock::now();

  auto dt_gen = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

  std::cout << "Generate: " << dt_gen.count() << " microseconds" << std::endl;

  t0 = std::chrono::high_resolution_clock::now();
  dispersion_filter(image, ref);
  t1 = std::chrono::high_resolution_clock::now();

  auto dt_ref = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

  std::cout << "Reference: " << dt_ref.count() << " microseconds" << std::endl;

  for (uint32_t i = 0, k = 0; i < NY; i++) {
    for (uint32_t j = 0; j < NX; j++, k++) {
      image_in[k] = image[k];
    }
  }

  auto t00 = std::chrono::high_resolution_clock::now();

  bo_in.sync(XCL_BO_SYNC_BO_TO_DEVICE);

  auto t01 = std::chrono::high_resolution_clock::now();

  auto run = krnl(bo_in, bo_out, NY);
  run.wait();

  auto t02 = std::chrono::high_resolution_clock::now();

  bo_out.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

  auto t03 = std::chrono::high_resolution_clock::now();

  auto t_in =
      std::chrono::duration_cast<std::chrono::microseconds>(t01 - t00).count();
  auto t_knl =
      std::chrono::duration_cast<std::chrono::microseconds>(t02 - t01).count();
  auto t_out =
      std::chrono::duration_cast<std::chrono::microseconds>(t03 - t02).count();

  std::cout << "FPGA: in / kernel / out: " << t_in << " / " << t_knl << " / "
            << t_out << " / "
            << " microseconds" << std::endl;

  for (int i = 0, k = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++, k++) {
      if (ref[k] != mask[k]) {
        wrong++;
        std::cout << "Pixel: " << k << " incorrect: " << ref[k] << " != " << mask[k] << std::endl;
      }
    }
  }

  if (wrong) {
    std::cout << wrong << " / " << ref.capacity() << " incorrect results"
              << std::endl;
  } else {
    std::cout << "Validation passed" << std::endl;
  }

  return (wrong > 0);
}
