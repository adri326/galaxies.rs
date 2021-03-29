# Galaxies.rs

A simple, multi-threaded, N-body simulation of galaxies.

**Demo:** https://www.youtube.com/watch?v=6vQRqB2M2yU

To speed up the process, this simulation uses an algorithm similar to the Barnes-Hut tree algorithm, by partitionning the space into two levels of voxel grids.

## Installation

```sh
git clone https://github.com/adri326/galaxies.rs
cd galaxies.rs
mkdir out/
# You may tweak some of the settings in src/main.rs, then run:
cargo run --release
# The window must currently be closed by hitting ctrl-c on the terminal
```

## TODO

- Argument to specify the output directory and automatically `mkdir` it (right now it has to be manually created)
- Actual Barnes-Hut tree algorithm (the current one has results close to it, but changing the scale of the simulation requires re-tweaking the constants to get the speed back)
- Pretty and efficient rendering (currently using `raqote` for rendering; all of the rendering happens in `src/draw.rs`)
- Add basic controls to the GUI (look around, close the window, etc.; right now the program has to be stopped with ctrl-c)
- Save/load states
- Quickly generate stable-ish galaxies (requires saving the state)
- Document and clean up the code (sorry~)
