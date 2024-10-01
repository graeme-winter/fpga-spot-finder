# FPGA Spot Finder

Spot finder, starting with the simple dispersion, built to demonstrate the code developed using the `vitis_hls` workflow [here](../019-project-fpga-summed-area/). This currently implements a streaming version of the spot finder compiled for a detector image the same size as an Eiger 2XE 16M: this is a compile time not run time option, since the hardware depends very significantly on the width of the image.

## Dependencies

The following are required to get anywhere at all:

- install of the Xilinx FPGA software tools (tested on 2022.2 versions)
- installation of XRT
- platform for the device you have which must support DMA (tested on U250 which is overkill)

The platform XPFM file must be passed either in the environment or via `make PLATFORM=`.

## Building

Building this is a multi-part process: build the kernel (`cd kernel; make PLATFORM=<platform>`), then link the kernel (`cd link; make PLATFORM=<platform>`), finally compile the host. Be warned that link the kernel takes hours, a chunk of CPU and a lot of memory. In all cases there is a `Makefile` which defines all the pertinent details, and all should _just work_.

## Code Design

The spot finder algorithm itself consists of three separate stages to compute the kernel functions on a stream rather than through local evaluation as performed on GPU. The stream components are:

- construct a summed area table for each of `m`, `i`, `i^2` - outputs a stream of blocks where each block is as wide as the image
- compute the local sum from the SAT for each of the parameters, over a square kernel region
- compute the threshold function based on the comparison of the pixel value with the local mean and variance

Around this there is a pixel doubler stage, to pass the pixel values into the summed area table calculation but also on to the thresholding with a suitable FIFO depth to allow the rest of the calculation to catch up. Finally, to map to the memory buffers, a reader and writer though I am not currently happy with these.

### Local Sum Evaluation

The summed area tables are typed as `uint8_t` for mask, `uint32_t` for sum of `i` and `i^2` - the latter of these is strictly incorrect for `uint16_t` data however:

- by the time the SATs are calculated the masked value pixels are set to zero
- the limiting factor is the sum over a square kernel region (i.e. 7x7 pixels)
- most data is generally smaller valued i.e. pixels > 1000 are _rare_

The SAT calculations could be split into three separate kernels, but they are currently implemented as a single kernel: the fact that all depend on computing the mask made it seem sensible however I am sure the conditional assignment on the pixel value is a minute amount of area. The input to the SAT accumulation is a stream of pixels, the output a stream of blocks - this made the local area calculation much simpler.

The summed area calculation is there to allow a local sum to be evaluated from four lookups:

```
D        C
 +-------+
 |       |
 |       |
 |   X   |
 |       |
 |       |
B+-------A
```

`A`, `B`, `C` and `D` are the sums from the upper left corner to that element in the array i.e. sum above and to the left of all the values in the array. The sum for the kernel region around `X` can be evaluated as `A + D - B - C`: this is simple and straightforward. And the most expensive part of the whole thing: this requires reading from four different positions in the array at the same time, while _writing_ to the same array in the same cycle with values coming in from the SAT calculation.

This was where the `hls::stream_of_blocks` became a lifesaver, because the `A` and `B` values can be read from the incoming block rather than this array, and written into this array from the block at the same time, ensuring an `II=1` can be maintained. This does however require a lot of area to deliver, and some partitioning of the SAT arrays, discussed below under "things I learned."

### Kernel Evaluation

Worth noting here that the filter calculation was recast to avoid any use of floating point values: everything is now done with integer arithmetic and no `sqrt` calls. The before and after:

```c++
      if (p > 0 && m_sum >= 2) {
        float bg_lhs = (float)m_sum * i2_sum - (float)i_sum * i_sum -
                       (float)i_sum * (m_sum - 1);
        float bg_rhs = i_sum * sigma_b * sqrtf((float)2.0f * (m_sum - 1));
        uint16_t background = bg_lhs > bg_rhs;
        float fg_lhs = (float)m_sum * p - (float)i_sum;
        float fg_rhs = sigma_s * sqrtf((float)i_sum * m_sum);
        uint16_t foreground = fg_lhs > fg_rhs;
        signal = background && foreground;
      }
```

```c++
        if (p > 0 && m_sum >= 2) {
          int64_t bg_lhs = (int64_t)m_sum * i2_sum - (int64_t)i_sum * i_sum -
                           (int64_t)i_sum * (m_sum - 1);
          int64_t bg_rhs2 = i_sum * i_sum * 72 * (m_sum - 1);
          uint8_t background = (bg_lhs > 0) && ((bg_lhs * bg_lhs) > bg_rhs2);
          int32_t fg_lhs = (int32_t)m_sum * p - (int32_t)i_sum;
          int32_t fg_rhs2 = 9 * i_sum * m_sum;
          uint8_t foreground = (fg_lhs > 0) && ((fg_lhs * fg_lhs) > fg_rhs2);
          signal = background && foreground;
        }
```

N.B. this compiles in `sigma_b = 6` and `sigma_s = 3` into the constant terms, and requires `uint64_t` for working variables for `uint16_t` measurements. The motivation behind this is the implementation of a floating point `hls::sqrtf` implies a latency of ~ 30 cycles, so to keep `N` pixels in flight and maintain `II=1` requires around `30*N` implementations of the calculation to keep the pipeline from stalling: a lot of area. As a aside, the term `sqrtf((float)2.0f * (m_sum - 1))` _can_ be tabulated since `m_sum - 1` has only around 50 possible valid values. The `bg_rhs` term can therefore be collapsed down to a lookup and a multiply.

## Debugging

Debugging the code turned out to be something of a nightmare in no small part because the behaviour of the system on the hardware is subtly but significantly different to the behaviour in the simulation flow in `vitis_hls`. To address this there are two additional steps that can be performed between the `vitis_hls` process and the hardware implementation: software and hardware emulation. The first of this is useful to check for "obvious" bugs (it is comparable to the `vitis_hls` process but includes simulation of the stream etc.) but this _does not_ really test everything. But it is quick.

### Software Emulation

In both the `kernel` and `link` directories, `make TARGET=sw_emu PLATFORM=<platform>` - this is relatively quick to execute. Then to run the emulation you need to pull a configuration from your platform with:

```
emconfigutil --platform <platform>
```

(this only needs to be done once to make the file) then run the simulation with:

```
XCL_EMULATION_MODE=sw_emu ../host/spot_finder ./spot_finder.xclbin
```

Yes, you need to have compiled the host for this to work. This runs pretty quickly (well, the kernel calculation is slow in comparison to the real kernel, but fast compared to `hw_emu` below) and should show that everything works as well as reporting the maximum depth of the simulation.

### Hardware Emulation

By analogy with the above, `make TARGET=hw_emu PLATFORM=<platform>` - again this is relatively quick to _compile_ compared with the real `hw` build. You also need the `emconfigutil` step if you didn't run this already, and may need a spot of `make clean` to remove the old `sw_emu` build targets.

```
XCL_EMULATION_MODE=hw_emu ../host/spot_finder ./spot_finder.xclbin
```

This will take a _long_ time to run - for example around 40 minutes to emulate 20ms of FPGA execution - but does so in a very verbose manner which is very useful for spotting bugs (e.g. deadlocks). There's a note on this under "things I learned."

## Things I learned

A lot of things were learned in the process, multiple days lost at various points trying to learn my way through the quagmire. Zeroth amongst these is the process embodied by `vitis_hls` is necessary but not sufficient for developing algorithms - it gives you some very useful guidance on how to write the code but guarantees nothing in terms of it actually working at the end of the day on real hardware.

### FPFA Development Boards

Developing only on huge data centre cards is _very slow_ because the compilation time takes hours. Embedded systems like Zynq 7000 are a lot smaller which makes for far faster compilation times - for testing ideas and learning your way around this is super helpful. The scripts here will work for smaller devices, provided you have a suitable platform file. Building the host will require a little more work, as will building a system image to deploy to a microSD card, which is slightly beyond the scope of this write up.

### Memory Usage

Algorithms need a lot more memory than you would think: at times (e.g. when building the local region calculation in stage 1) a lot of reads happen concurrently from a shared buffer: the compilers make this possilble by maintaining _multiple copies_ of that buffer which all occupy BRAM or area. This means, for something where you would expect to use e.g. 14 x 18 kb BRAMs you may find it using 4 or 8 times as many to allow the read concurrency in addition to specified partitioning. This is at the root of why the calculation for one Eiger 2XE module (512 rows of 1028 pixels) won't fit onto the Zynq 7020 in a Zybo Z7-20 board.

### Stream Naming

This one was a real doozy and only exposed itself when trying to run on real hardware (or in `hw_emu` mode). Spoilers: HLS may _look_ like C++ and use C++ syntax but it _does not_ follow all C++ rules. For example in:

```c++
void function(type_stream_t &in, type_stream_t &out) {
  // implementation
}
```

`in` and `out` are _in scope_ i.e. these are locally defined and any reference to these is independent of the same names in other contexts. THIS IS NOT TRUE IN HLS. However nothing in the software chain _tells_ you that this is not true, and all the software emulation parts above in `vitis_hls` work ✨just fine✨.

Until you hit a deadlock and you have no idea why. This was where the `make TARGET=hw_emu` etc. above was überhelful because it _actually crashed_. When you run in the test environment a lot of debug information is output to a text file `output.txt` in a directory reported at the start of the output. Searching for `ERROR` in here showed:

```
//////////////////////////////////////////////////////////////////////////////
// ERROR!!! DEADLOCK DETECTED at 217790000 ns! SIMULATION WILL BE STOPPED! //
//////////////////////////////////////////////////////////////////////////////
/////////////////////////
// Dependence cycle 1:
// (1): Process: spot_finder_spot_finder.spot_finder_kernel_U0.duplicate_pack16_stream_U0
//      Blocked by empty input FIFO 'spot_finder_spot_finder.input_U' written by process 'spot_finder_spot_finder.reader_U0'
// (2): Process: spot_finder_spot_finder.reader_U0
//      Blocked by empty input FIFO 'spot_finder_spot_finder.input_U' written by process 'spot_finder_spot_finder.spot_finder_kernel_U0.duplicate_pack16_stream_U0'
////////////////////////////////////////////////////////////////////////
// Totally 1 cycles detected!
////////////////////////////////////////////////////////////////////////
// ERROR!!! DEADLOCK DETECTED at 217790000 ns! SIMULATION WILL BE STOPPED! //
//////////////////////////////////////////////////////////////////////////////
```

Wait, says I: the input and output to this were _the same_?

This is because I was being lazy and making up a pattern where the reader and writer both just called the input and output channels `input` and `output` - as well as other channels for e.g. the pixel splitting path. The actual implementation aliased these because there is no scope & hence broke. Hacking in some ugly names made the emulation actually run, then the real version actually ran.

## Results

Happily I survived the learning process and debugging and got to some working code. The test case shows the calculation being reasonably quick, but not really worth the effort at the moment:

```
cs04r-sc-serv-85 host :) [run-fpga-spot] $ ./spot_finder ../link/spot_finder.xclbin
Generate: 3899034 microseconds
Reference: 231416 microseconds
Done copying...
FPGA: in / kernel / out: 3663 / 28484 / 3097 /  microseconds
0 / 17919096 incorrect results
```

28 ms to perform the calculation equates to around 35 frames / second - a very long way short of the target 500 frames / second. Overall `II=4510889` cycles, so average `II=1` for 4 pixel blocks maintiained. In implementation the minimum cycle time was 5.34ns which explains the slow evaluation.

## Future Performance Improvement

The current implementation takes packs of four pixels and performs the corresponding calculations concurrently - so for an image of 4362 x 4108 pixels (an Eiger 2 XE 16M) we would expect a minimum of around 15 ms to execute at 300 MHz - the kernel execution time of 28 ms is therefore already excessive, in addition to being fundamentally slow.

Ingesting more pixels per clock (i.e. `PACK=8, 16, 32`) would immediately improve the performance however the row data could need to be packed with NULL values to keep the width a multiple of `PACK`. There will also be performance tuning implications of doing this, e.g. on the lengths of the most delayed paths through the calculation.

Running multiple _instances_ of the kernel, to analyse images concurrently, or partitioning the image into e.g. 32 modules and computing concurrently, would also be opportunities to improve the throughput. Performing the image decomopression on the card unlikely to be significant, as the transfer time is around 3-4 ms, but for multiple images pipelining the calculation to overlap the transfer with kernel execution, as performed for the GPU equivalent, would be worthwhile.
