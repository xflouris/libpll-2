  ocho $10
  ocho t are already using libpll-2
   List of projects already using libpll-2 and site repeats, and reported speedups compared with the tip pattern optimization:
   * [RAxML-NG](https://github.com/amkozlov/raxml-ng): speedup ranges between 1.2 and 1.5 
   * [ModelTest-NG](https://github.com/ddarriba/modeltest): speedup around 2
   * [EPA-ng](https://github.com/Pbdas/epa-ng): no speedup (likelihood computation is not the main bottleneck) but memory footprint reduced by 30%.$10
# Libpll-t are already using libpll-2
   List of projects already using libpll-2 and site repeats, and reported speedups compared with the tip pattern optimization:
   * [RAxML-NG](https://github.com/amkozlov/raxml-ng): speedup ranges between 1.2 and 1.5 
   * [ModelTest-NG](https://github.com/ddarriba/modeltest): speedup around 2
   * [EPA-ng](https://github.com/Pbdas/epa-ng): no speedup (likelihood computation is not the main bottleneck) but memory footprint reduced by 30%.2

libpll-2 is the new official fork of libpll (https://github.com/xflouris/libpll/). It implements site repeats to speed up computations.


Please read the wiki for more information.



# Projects that are already using libpll-2
 List of projects already using libpll-2 and site repeats, and reported speedups compared with the tip pattern optimization:
 * [RAxML-NG](https://github.com/amkozlov/raxml-ng): speedup ranges between 1.2 and 1.5 
 * [ModelTest-NG](https://github.com/ddarriba/modeltest): speedup around 2
 * [EPA-ng](https://github.com/Pbdas/epa-ng): no speedup (likelihood computation is not the main bottleneck) but memory footprint reduced by 30%.


# Compilation instructions

Currently, `libpll` requires that [GNU Bison](http://www.gnu.org/software/bison/)
and [Flex](http://flex.sourceforge.net/) are installed on the target system. On
a Debian-based Linux system, the two packages can be installed using the command

`apt-get install flex bison`

The library also requires that a GNU system is available as it uses several
functions (e.g. `asprintf`) which are not present in the POSIX standard.
This, however will change in the future in order to have a more portable
and cross-platform library.

The library can be compiled using either of the following two ways.

**Cloning the repo** Clone the repo and bild the executable and documentation
using the following commands.

```bash
git clone https://github.com/xflouris/libpll.git
cd libpll
./autogen.sh
./configure
make
make install    # as root, otherwise run: sudo make install
```

When using the cloned repository version, you will also need
[autoconf](https://www.gnu.org/software/autoconf/autoconf.html),
[automake](https://www.gnu.org/software/automake/) and
[libtool](https://www.gnu.org/software/libtool/) installed. On a Debian-based
Linux system, the packages can be installed using the command

```bash
sudo apt-get install autotools-dev autoconf libtool
```

The library will be installed on the operating system's standard paths.  For
some GNU/Linux distributions it might be necessary to add that standard path
(typically `/usr/local/lib`) to `/etc/ld.so.conf` and run `ldconfig`.

Microsoft Windows compatibility was tested with a cross-compiler and seems to
work out-of-the-box using [MingW](http://www.mingw.org/).

# libpll-2 license and third party licenses

The libpll-2 code is currently licensed under the
[GNU Affero General Public License version 3](http://www.gnu.org/licenses/agpl-3.0.en.html).
Please see LICENSE.txt for details.

libpll-2 includes code from several other projects. We would like to thank the
authors for making their source code available.

libpll includes code from GNU Compiler Collection distributed under the GNU
General Public License.



