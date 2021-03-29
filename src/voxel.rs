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

    pub children: Vec<Voxel>,
}

impl Voxel {
    pub fn from_pos(resolution: Prec, position: (Prec, Prec, Prec), index: usize, children: Vec<Voxel>) -> Self {
        Self {
            coords: get_coords(resolution, position),
            position: round_resolution(resolution, position),
            elements: vec![index],
            total_mass: 0.0,
            center_of_mass: (0.0, 0.0, 0.0),
            average_speed: (0.0, 0.0, 0.0),
            children,
        }
    }
}

#[derive(Debug, Clone)]
pub struct VoxelGrid {
    pub voxels: Vec<Voxel>,
    pub grid: HashMap<(i32, i32, i32), usize>,

    pub resolution: Prec,
    pub sub_resolution: Prec,

    pub gone_radius: Prec,
}

impl VoxelGrid {
    pub fn new(resolution: Prec, sub_resolution: Prec, gone_radius: Prec) -> Self {
        Self {
            voxels: Vec::new(),
            grid: HashMap::new(),

            resolution,
            sub_resolution,
            gone_radius,
        }
    }

    pub fn populate(
        &mut self,
        position: &Table<(Prec, Prec, Prec)>,
    ) {
        for (index, &p) in position.iter().enumerate() {
            let coords = get_coords(self.resolution, p);

            if distance_squared(p, (0.0, 0.0, 0.0)) > self.gone_radius * self.gone_radius {
                continue
            }

            if let Some(&i) = self.grid.get(&coords) {
                self.voxels[i].elements.push(index);
                // Append to the sub-voxels
                let mut found_sub_voxel = false;
                let sub_coords = get_coords(self.sub_resolution, p);

                for sub_voxel in self.voxels[i].children.iter_mut() {
                    if sub_voxel.coords == sub_coords {
                        sub_voxel.elements.push(index);
                        found_sub_voxel = true;
                    }
                }

                if !found_sub_voxel {
                    self.voxels[i].children.push(
                        Voxel::from_pos(self.sub_resolution, p, index, vec![])
                    );
                }
            } else {
                self.voxels.push(
                    Voxel::from_pos(self.resolution, p, index, vec![
                        Voxel::from_pos(self.sub_resolution, p, index, vec![]),
                    ])
                );
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

            if total_mass > 0.0 {
                voxel.center_of_mass = mul(position_sum, 1.0 / total_mass);
                voxel.total_mass = total_mass;
                voxel.average_speed = mul(speed_sum, 1.0 / total_mass);
            } else {
                voxel.center_of_mass = voxel.position;
                voxel.total_mass = total_mass;
                voxel.average_speed = (0.0, 0.0, 0.0);
            }

            for sub_voxel in &mut voxel.children {
                let mut position_sum = (0.0, 0.0, 0.0);
                let mut speed_sum = (0.0, 0.0, 0.0);
                let mut total_mass = 0.0;

                for &i in &sub_voxel.elements {
                    let p = *position.get(i).unwrap();
                    let s = *speed.get(i).unwrap();
                    let m = *mass.get(i).unwrap();
                    position_sum = add(position_sum, mul(p, m));
                    speed_sum = add(speed_sum, mul(s, m));
                    total_mass += m;
                }

                if total_mass > 0.0 {
                    sub_voxel.center_of_mass = mul(position_sum, 1.0 / total_mass);
                    sub_voxel.total_mass = total_mass;
                    sub_voxel.average_speed = mul(speed_sum, 1.0 / total_mass);
                } else {
                    sub_voxel.center_of_mass = sub_voxel.position;
                    sub_voxel.total_mass = 0.0;
                    sub_voxel.average_speed = (0.0, 0.0, 0.0);
                }
            }
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
                let sub_coords = get_coords(self.sub_resolution, *position.get(index).unwrap());

                if coords != voxel.coords {
                    voxel.elements.swap_remove(i);
                    to_relocate.push((index, coords));

                    // Remove from sub-voxel
                    for sub_voxel in voxel.children.iter_mut() {
                        if sub_voxel.coords == sub_coords {
                            let mut j = 0;
                            while j < sub_voxel.elements.len() {
                                if sub_voxel.elements[j] == index {
                                    sub_voxel.elements.swap_remove(j);
                                    break
                                } else {
                                    j += 1;
                                }
                            }
                        }
                    }
                } else {
                    // Relocate the star within the sub-voxels if needed
                    for sub_voxel in voxel.children.iter_mut() {
                        if sub_voxel.coords == sub_coords {
                            if !sub_voxel.elements.iter().any(|&e| e == index) {
                                sub_voxel.elements.push(index);
                            }
                        } else {
                            let mut j = 0;
                            while j < sub_voxel.elements.len() {
                                if sub_voxel.elements[j] == index {
                                    sub_voxel.elements.swap_remove(j);
                                    break
                                } else {
                                    j += 1;
                                }
                            }
                        }
                    }


                    i += 1;
                }
            }
        }

        for (index, coords) in to_relocate {
            let p = *position.get(index).unwrap();
            if distance_squared(p, (0.0, 0.0, 0.0)) > self.gone_radius * self.gone_radius {
                continue
            }
            if let Some(&i) = self.grid.get(&coords) {
                self.voxels[i].elements.push(index);
                // Append to the sub-voxels
                let mut found_sub_voxel = false;
                let sub_coords = get_coords(self.sub_resolution, p);

                for sub_voxel in self.voxels[i].children.iter_mut() {
                    if sub_voxel.coords == sub_coords {
                        sub_voxel.elements.push(index);
                        found_sub_voxel = true;
                    }
                }

                if !found_sub_voxel {
                    self.voxels[i].children.push(
                        Voxel::from_pos(self.sub_resolution, p, index, vec![])
                    );
                }
            } else {
                self.voxels.push(Voxel::from_pos(self.resolution, p, index, vec![
                    Voxel::from_pos(self.sub_resolution, p, index, vec![]),
                ]));
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
                let mut j = 0;
                while j < self.voxels[i].children.len() {
                    if self.voxels[i].children[j].elements.len() == 0 {
                        self.voxels[i].children.swap_remove(j);
                    } else {
                        j += 1;
                    }
                }
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

pub fn round_resolution(resolution: Prec, position: (Prec, Prec, Prec)) -> (Prec, Prec, Prec) {
    (
        (position.0 / resolution).round() * resolution,
        (position.1 / resolution).round() * resolution,
        (position.2 / resolution).round() * resolution,
    )
}
