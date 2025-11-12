# libAnnalee

Marko Grönroos 1998-2005.

libAnnalee is a neuroevolution library for research purposes.
It uses an object-oriented neural model that allows complex structures.

The library requires [MagiCLib++](/magi42/magiclib) common base library,
[libinanna](/magi42/libinanna) Articial Neural Network library, and
[libnhp](/magi42/libinanna) evolutionary algorithm library.

It is expected to be compiled under the MagiCLib++ source tree, to be able to use and develop the base library more easily.
To have it compile there, you need to include it in the `Makefile` of MagiCLib++, with the `makemodules` parameter:

```
makemodules = libmagic libinanna libnhp libannalee
```

Then, you can run `make` to compile the dependencies and the library itself.

## Projects

Under the `projects`, there are applications using the library.
