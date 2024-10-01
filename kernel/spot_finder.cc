#include <hls_math.h>

#include "spot_finder.h"

void duplicate_pack16_stream(pack16_stream &source, pack16_stream &output0,
                             pack16_stream &output1, int rows) {
#pragma HLS STREAM type = fifo depth = (NX * 2) variable = output1
  for (int32_t i = 0; i < rows; i++) {
    pack16_t p;
    for (int j = 0; j < NX; j += PACK) {
#pragma HLS PIPELINE II = 1
      p = source.read();
      output0.write(p);
      output1.write(p);
    }
  }
}

void stage0(pack16_stream &input, block8_stream &m_sat_stream,
            block32_stream &i_sat_stream, block32_stream &i2_sat_stream,
            int rows) {
  uint8_t acc_m[NX] = {0};
  uint32_t acc_i[NX] = {0};
  uint32_t acc_i2[NX] = {0};
#pragma HLS ARRAY_PARTITION variable = acc_m type = cyclic factor = PACK
#pragma HLS ARRAY_PARTITION variable = acc_i type = cyclic factor = PACK
#pragma HLS ARRAY_PARTITION variable = acc_i2 type = cyclic factor = PACK

sat_row_loop:
  for (int32_t row = 0; row < rows; row++) {
    hls::write_lock<block8_t> m_lock(m_sat_stream);
    hls::write_lock<block32_t> i_lock(i_sat_stream);
    hls::write_lock<block32_t> i2_lock(i2_sat_stream);
    uint8_t line_m = 0;
    uint32_t line_i = 0, line_i2 = 0;
    pack16_t p;
    pack8_t m, m_out;
    pack32_t i, i2, i_out, i2_out;
  sat_pack_loop:
    for (uint32_t pxl = 0; pxl < NX; pxl += PACK) {
#pragma HLS PIPELINE II = 1
      p = input.read();
      uint8_t scr_m[PACK];
      uint32_t scr_i[PACK];
      uint32_t scr_i2[PACK];
#pragma HLS ARRAY_PARTITION variable = scr_m type = complete
#pragma HLS ARRAY_PARTITION variable = scr_i type = complete
#pragma HLS ARRAY_PARTITION variable = scr_i2 type = complete
    sat_pxl_loop_a:
      for (int l0 = 0; l0 < PACK; l0++) {
#pragma HLS UNROLL
        m[l0] = p[l0] > 0xfffd ? 0 : 1;
        i[l0] = p[l0] > 0xfffd ? 0 : p[l0];
        i2[l0] = p[l0] > 0xfffd ? 0 : p[l0] * p[l0];
      }
    sat_pxl_loop_b:
      for (int l0 = 0; l0 < PACK; l0++) {
#pragma HLS UNROLL
        scr_m[l0] = 0;
        scr_i[l0] = 0;
        scr_i2[l0] = 0;
      sat_pxl_loop_c:
        for (int l1 = 0; l1 <= l0; l1++) {
#pragma HLS UNROLL
          scr_m[l0] += m[l1];
          scr_i[l0] += i[l1];
          scr_i2[l0] += i2[l1];
        }
        acc_m[pxl + l0] += line_m + scr_m[l0];
        acc_i[pxl + l0] += line_i + scr_i[l0];
        acc_i2[pxl + l0] += line_i2 + scr_i2[l0];
      }
      for (int l0 = 0; l0 < PACK; l0++) {
#pragma HLS UNROLL
        m_out[l0] = acc_m[pxl + l0];
        i_out[l0] = acc_i[pxl + l0];
        i2_out[l0] = acc_i2[pxl + l0];
      }
      m_lock[pxl / PACK] = m_out;
      i_lock[pxl / PACK] = i_out;
      i2_lock[pxl / PACK] = i2_out;

      line_m += scr_m[PACK - 1];
      line_i += scr_i[PACK - 1];
      line_i2 += scr_i2[PACK - 1];
    }
  }
}

void stage1(block8_stream &m_sat_stream, block32_stream &i_sat_stream,
            block32_stream &i2_sat_stream, block8_stream &m_loc_stream,
            block32_stream &i_loc_stream, block32_stream &i2_loc_stream,
            int rows) {

  uint8_t t_m[KNL2][NX];
  uint32_t t_i[KNL2][NX];
  uint32_t t_i2[KNL2][NX];
#pragma HLS ARRAY_PARTITION variable = t_m type = cyclic dim = 2 factor = PACK
#pragma HLS ARRAY_PARTITION variable = t_i type = cyclic dim = 2 factor = PACK
#pragma HLS ARRAY_PARTITION variable = t_i2 type = cyclic dim = 2 factor = PACK

pre_row_loop:
  for (int32_t row = 0; row < KNL; row++) {
    hls::read_lock<block8_t> m_slock(m_sat_stream);
    hls::read_lock<block32_t> i_slock(i_sat_stream);
    hls::read_lock<block32_t> i2_slock(i2_sat_stream);

    uint32_t t_row = row;

  pre_pack_loop:
    for (uint32_t pxl = 0; pxl < NX; pxl += PACK) {
#pragma HLS PIPELINE II = 1
      pack8_t m, m_out;
      pack32_t i, i2, i_out, i2_out;

      m = m_slock[pxl / PACK];
      i = i_slock[pxl / PACK];
      i2 = i2_slock[pxl / PACK];
    pre_pxl_loop:
      for (int l0 = 0; l0 < PACK; l0++) {
#pragma HLS UNROLL
        t_m[t_row][pxl + l0] = m[l0];
        t_i[t_row][pxl + l0] = i[l0];
        t_i2[t_row][pxl + l0] = i2[l0];
      }
    }
  }

diff_row_loop:
  for (int32_t row = KNL; row < rows; row++) {
    hls::read_lock<block8_t> m_slock(m_sat_stream);
    hls::read_lock<block32_t> i_slock(i_sat_stream);
    hls::read_lock<block32_t> i2_slock(i2_sat_stream);
    hls::write_lock<block8_t> m_llock(m_loc_stream);
    hls::write_lock<block32_t> i_llock(i_loc_stream);
    hls::write_lock<block32_t> i2_llock(i2_loc_stream);
    int32_t t_row = row % KNL2;

    int32_t i0 = row - 2 * KNL - 1;
    int32_t i1 = row;

    int32_t ti0 = i0 % KNL2;
    int32_t ti1 = i1 % KNL2;

  diff_pack_loop:
    for (uint32_t pxl = 0; pxl < NX; pxl += PACK) {
#pragma HLS PIPELINE II = 1

      pack8_t m;
      pack32_t i, i2;

      pack8_t m_out;
      pack32_t i_out, i2_out;

      m = m_slock[pxl / PACK];
      i = i_slock[pxl / PACK];
      i2 = i2_slock[pxl / PACK];

    diff_pxl_loop:
      for (int l0 = 0; l0 < PACK; l0++) {
#pragma HLS UNROLL
        int j = pxl + l0;
        int j0 = j - KNL - 1;
        int j1 = j + KNL >= NX ? NX - 1 : j + KNL;

        t_m[ti1][j] = m[l0];
        t_i[ti1][j] = i[l0];
        t_i2[ti1][j] = i2[l0];

        uint8_t m_sum_a = m_slock[j1 / PACK][j1 & (PACK - 1)];
        uint32_t i_sum_a = i_slock[j1 / PACK][j1 & (PACK - 1)];
        uint32_t i2_sum_a = i2_slock[j1 / PACK][j1 & (PACK - 1)];

        uint8_t m_sum_b = j0 >= 0 ? m_slock[j0 / PACK][j0 & (PACK - 1)] : 0;
        uint32_t i_sum_b = j0 >= 0 ? i_slock[j0 / PACK][j0 & (PACK - 1)] : 0;
        uint32_t i2_sum_b = j0 >= 0 ? i2_slock[j0 / PACK][j0 & (PACK - 1)] : 0;

        uint8_t m_sum_c = i0 >= 0 ? t_m[ti0][j1] : 0;
        uint32_t i_sum_c = i0 >= 0 ? t_i[ti0][j1] : 0;
        uint32_t i2_sum_c = i0 >= 0 ? t_i2[ti0][j1] : 0;

        uint8_t m_sum_d = (j0 >= 0 && i0 >= 0) ? t_m[ti0][j0] : 0;
        uint32_t i_sum_d = (j0 >= 0 && i0 >= 0) ? t_i[ti0][j0] : 0;
        uint32_t i2_sum_d = (j0 >= 0 && i0 >= 0) ? t_i2[ti0][j0] : 0;

        m_out[l0] = m_sum_a + m_sum_d - m_sum_b - m_sum_c;
        i_out[l0] = i_sum_a + i_sum_d - i_sum_b - i_sum_c;
        i2_out[l0] = i2_sum_a + i2_sum_d - i2_sum_b - i2_sum_c;
      }

      m_llock[pxl / PACK] = m_out;
      i_llock[pxl / PACK] = i_out;
      i2_llock[pxl / PACK] = i2_out;
    }
  }

post_row_loop:
  for (int32_t row = rows; row < rows + KNL; row++) {
    hls::write_lock<block8_t> m_llock(m_loc_stream);
    hls::write_lock<block32_t> i_llock(i_loc_stream);
    hls::write_lock<block32_t> i2_llock(i2_loc_stream);
    int32_t t_row = row % KNL2;

    int32_t i0 = row - 2 * KNL - 1;
    int32_t i1 = rows + KNL;

    int32_t ti0 = i0 % KNL2;
    int32_t ti1 = i1 % KNL2;

  post_pack_loop:
    for (uint32_t pxl = 0; pxl < NX; pxl += PACK) {
#pragma HLS PIPELINE II = 1
      pack8_t m, m_out;
      pack32_t i, i2, i_out, i2_out;
    post_pxl_loop:
      for (int l0 = 0; l0 < PACK; l0++) {
#pragma HLS UNROLL
        int j = pxl + l0;
        int j0 = j - KNL - 1;
        int j1 = j + KNL >= NX ? NX - 1 : j + KNL;

        uint8_t m_sum_ac = t_m[ti1][j1] - t_m[ti0][j1];
        uint32_t i_sum_ac = t_i[ti1][j1] - t_i[ti0][j1];
        uint32_t i2_sum_ac = t_i2[ti1][j1] - t_i2[ti0][j1];

        uint8_t m_sum_bd = j0 >= 0 ? t_m[ti0][j0] - t_m[ti1][j0] : 0;
        uint32_t i_sum_bd = j0 >= 0 ? t_i[ti0][j0] - t_i[ti1][j0] : 0;
        uint32_t i2_sum_bd = j0 >= 0 ? t_i2[ti0][j0] - t_i2[ti1][j0] : 0;

        m_out[l0] = m_sum_ac + m_sum_bd;
        i_out[l0] = i_sum_ac + i_sum_bd;
        i2_out[l0] = i2_sum_ac + i2_sum_bd;
      }

      m_llock[pxl / PACK] = m_out;
      i_llock[pxl / PACK] = i_out;
      i2_llock[pxl / PACK] = i2_out;
    }
  }
}

void stage2(pack16_stream &p_stream, block8_stream &m_loc_stream,
            block32_stream &i_loc_stream, block32_stream &i2_loc_stream,
            pack16_stream &output, int rows) {

filter_row_loop:
  for (int32_t row = 0; row < rows; row++) {
    hls::read_lock<block8_t> m_lock(m_loc_stream);
    hls::read_lock<block32_t> i_lock(i_loc_stream);
    hls::read_lock<block32_t> i2_lock(i2_loc_stream);
    pack8_t m;
    pack16_t pack, q;
    pack32_t i, i2;
  filter_pack_loop:
    for (uint32_t pxl = 0; pxl < NX; pxl += PACK) {
#pragma HLS PIPELINE II = 1
      pack = p_stream.read();
      m = m_lock[pxl / PACK];
      i = i_lock[pxl / PACK];
      i2 = i2_lock[pxl / PACK];
    filter_pxl_loop:
      for (int l0 = 0; l0 < PACK; l0++) {
#pragma HLS UNROLL
        uint8_t m_sum = m[l0];
        uint32_t i_sum = i[l0];
        uint32_t i2_sum = i2[l0];
        uint16_t p = pack[l0] > 0xfffd ? 0 : pack[l0];

        uint8_t s = 0;

        if (p > 0 && m_sum >= 2) {
          int64_t bg_lhs = (int64_t)m_sum * i2_sum - (int64_t)i_sum * i_sum -
                           (int64_t)i_sum * (m_sum - 1);
          int64_t bg_rhs2 = i_sum * i_sum * 72 * (m_sum - 1);
          uint8_t background = (bg_lhs > 0) && ((bg_lhs * bg_lhs) > bg_rhs2);
          int32_t fg_lhs = (int32_t)m_sum * p - (int32_t)i_sum;
          int32_t fg_rhs2 = 9 * i_sum * m_sum;
          uint8_t foreground = (fg_lhs > 0) && ((fg_lhs * fg_lhs) > fg_rhs2);
          s = background && foreground;
        }

        q[l0] = s ? p : 0;
      }
      output.write(q);
    }
  }
}

void reader(pack16_t *read_input, pack16_stream &read_output, int rows) {
  for (int row = 0, k = 0; row < rows; row++) {
    for (int pack = 0; pack < NX; pack += PACK, k++) {
      read_output.write(read_input[k]);
    }
  }
}

void writer(pack16_stream &write_input, pack16_t *write_output, int rows) {
  for (int row = 0, k = 0; row < rows; row++) {
    for (int pack = 0; pack < NX; pack += PACK, k++) {
      write_output[k] = write_input.read();
    }
  }
}

void spot_finder_kernel(pack16_stream &kernel_input, pack16_stream &kernel_output, int rows) {
#pragma HLS DATAFLOW
  pack16_stream p_stream0, p_stream1;
  block8_stream m_sat_stream, m_loc_stream;
  block32_stream i_sat_stream, i_loc_stream, i2_sat_stream, i2_loc_stream;

  duplicate_pack16_stream(kernel_input, p_stream0, p_stream1, rows);

  stage0(p_stream0, m_sat_stream, i_sat_stream, i2_sat_stream, rows);
  stage1(m_sat_stream, i_sat_stream, i2_sat_stream, m_loc_stream, i_loc_stream,
         i2_loc_stream, rows);
  stage2(p_stream1, m_loc_stream, i_loc_stream, i2_loc_stream, kernel_output, rows);
}

extern "C" {
void spot_finder(pack16_t *top_in, pack16_t *top_out, int rows) {
  static pack16_stream top_input;
  static pack16_stream top_output;

#pragma HLS INTERFACE m_axi port = top_in bundle = gmem0
#pragma HLS INTERFACE m_axi port = top_out bundle = gmem1
#pragma HLS DATAFLOW

  reader(top_in, top_input, rows);
  spot_finder_kernel(top_input, top_output, rows);
  writer(top_output, top_out, rows);
}
}