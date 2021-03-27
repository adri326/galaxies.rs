use super::*;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Voxel {
    pub coords: (i32, i32, i32),
    pub position: (Prec, Prec, Prec),
    // pub neihbors: [Option<usize>; 26],
    pub elements: Vec<usize>,
    pub total_mass: Prec,
    pub center_of_mass: (Prec, Prec, Prec),
    pub average_speed: (Prec, Prec, Prec),
}

pub struct VoxelGrid {
    pub voxels: Vec<Voxel>,
    pub grid: HashMap<(i32, i32, i32), usize>,

    pub resolution: Prec,
}

impl VoxelGrid {
    pub fn new(resolution: Prec) -> Self {
        Self {
            voxels: Vec::new(),
            grid: HashMap::new(),

            resolution,
        }
    }

    pub fn populate(
        &mut self,
        position: &Table<(Prec, Prec, Prec)>,
    ) {
        for (index, &p) in position.iter().enumerate() {
            let coords = get_coords(self.resolution, p);
            if let Some(&i) = self.grid.get(&coords) {
                self.voxels[i].elements.push(index);
            } else {
                self.voxels.push(Voxel {
                    coords,
                    position: get_position(self.resolution, coords),
                    elements: vec![index],
                    total_mass: 0.0,
                    center_of_mass: (0.0, 0.0, 0.0),
                    average_speed: (0.0, 0.0, 0.0),
                });
                self.grid.insert(coords, self.voxels.len() - 1);
            }
        }
    }

    pub fn update(
        &mut self,
        position: &Table<(Prec, Prec, Prec)>,
        speed: &Table<(Prec, Prec, Prec)>,
        mass: &Table<Prec>,
    ) {
        for voxel in &mut self.voxels {
            let mut position_sum = (0.0, 0.0, 0.0);
            let mut speed_sum = (0.0, 0.0, 0.0);
            let mut total_mass = 0.0;

            for &i in &voxel.elements {
                let p = *position.get(i).unwrap();
                let s = *speed.get(i).unwrap();
                let m = *mass.get(i).unwrap();
                position_sum = add(position_sum, mul(p, m));
                speed_sum = add(speed_sum, mul(s, m));
                total_mass += m;
            }

            voxel.center_of_mass = mul(position_sum, 1.0 / total_mass);
            voxel.total_mass = total_mass;
            voxel.average_speed = mul(speed_sum, 1.0 / total_mass);
        }
    }

    pub fn relocate(
        &mut self,
        position: &Table<(Prec, Prec, Prec)>,
    ) {
        let mut to_relocate = Vec::new();
        for voxel in &mut self.voxels {
            let mut i = 0;
            while i < voxel.elements.len() {
                let index = voxel.elements[i];
                let coords = get_coords(self.resolution, *position.get(index).unwrap());
                if coords != voxel.coords {
                    voxel.elements.swap_remove(i);
                    to_relocate.push((index, coords));
                } else {
                    i += 1;
                }
            }
        }

        for (index, coords) in to_relocate {
            if let Some(&i) = self.grid.get(&coords) {
                self.voxels[i].elements.push(index);
            } else {
                self.voxels.push(Voxel {
                    coords,
                    position: get_position(self.resolution, coords),
                    elements: vec![index],
                    total_mass: 0.0,
                    center_of_mass: (0.0, 0.0, 0.0),
                    average_speed: (0.0, 0.0, 0.0),
                });
                self.grid.insert(coords, self.voxels.len() - 1);
            }
        }
    }

    pub fn gc(&mut self) {
        let mut i = 0;
        while i < self.voxels.len() {
            if self.voxels[i].elements.len() == 0 {
                self.grid.remove(&self.voxels[i].coords);
                self.voxels.swap_remove(i);
            } else {
                i += 1;
            }
        }

        for (index, voxel) in self.voxels.iter().enumerate() {
            self.grid.insert(voxel.coords, index);
        }
    }

    pub fn iter<'a>(&'a self) -> std::slice::Iter<'a, Voxel> {
        self.voxels.iter()
    }
}

pub fn get_coords(resolution: Prec, position: (Prec, Prec, Prec)) -> (i32, i32, i32) {
    (
        (position.0 / resolution).round() as i32,
        (position.1 / resolution).round() as i32,
        (position.2 / resolution).round() as i32,
    )
}

pub fn get_position(resolution: Prec, coords: (i32, i32, i32)) -> (Prec, Prec, Prec) {
    (
        (coords.0 as Prec * resolution),
        (coords.1 as Prec * resolution),
        (coords.2 as Prec * resolution),
    )
}
