# x87 library

This library is designed to give authors of x86 emulators a resource that can be used for implementing x87 operations.
It is currently used for x87 floating-point support in my DREAMM emulator.

The motivation behind this is that existing libraries providing software floating point operations tend to be general-purpose and require a lot of tweaking/massaging in order to make them viable to use for x87 emulation.
The goal of this library is provide a set of tested implementations that not only produce accurate results but also generate status flags and handle edge cases in ways that match x87 implementations.

The core headers are **x87fp64.h** and **x87fp80.h**, which provide two types **x87::fp64_t** and **x87::fp80_t** respectively.
These two types have identical interfaces and are intended to be easily swappable depending on your needs.
**x87::fp64_t** performs all math using native 64-bit double support on the current processor (assumes either x64 or ARM64).
**x87::fp80_t** by contrast performs all math operations by hand to full 80-bit precision.

Please note, however, that at this time the full 80-bit code has not been implemented, so really **x87::fp64_t** is the only complete implementation. **x87::fp80_t** does, however, have a well-tested set of loads and stores, including integer conversions.

Feel free to use this code in your projects if it is useful.
And if you find any bugs or the motivation to enhance/improve it in any way, I am definitely open to any improvements!
