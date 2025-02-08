#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]

use minifb::{Window, WindowOptions};
use rand::Rng;
use std::fs::File;
use std::io::{self, Write};
use std::{thread, time};

use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

//Simulation settings
const NUM_RODS: usize = 20;
const NUM_BOXES: usize = 5;
const MAX_BOXES: usize = 12; //you can only have so many
const MAX_RODS: usize = GRID_SIZE*GRID_SIZE / (2*2);
const GRID_SIZE: usize = 30; //2 * (MOON_RADIUS as usize);
const NUM_GENERATIONS: u32 = 1000;
const POPULATION_SIZE: u64 = 100;

//Graphics settings
const PARTICLE_RADIUS: usize = SCALE / 8;
const WIDTH: usize = 1000; //GRID_SIZE * SCALE;
const HEIGHT: usize = WIDTH; //GRID_SIZE * SCALE;
const SCALE: usize = ((WIDTH as f64) / (GRID_SIZE as f64)) as usize;
const FRAME_TIME: time::Duration = time::Duration::from_millis(16); // ~60 FPS

// Simulation constants - Known
// const MOON_RADIUS: f64 = 142.0;
const SOLAR_WIND_MIN_VELOCITY: f64 = 0.225;

const ROD_DEFLECTION_STRENGTH: f64 = 4.0;
const ROD_MAX_RANGE_SQUARED: f64 = 25.0 * 25.0;
const MIN_ELECTROMAGNETIC_DISTANCE: f64 = 2.0;
const SOLAR_WIND_DEFLECTION_TICK_INTERVAL: u32 = 6;
const SOLAR_WIND_TICK_MULTIPLIER: u32 = 1; //1 when watching, 12 when not?
const COOLDOWN_DISTANCE: f64 = 1.5;
const COOLDOWN_TICKS: f64 = 30.0;

// Changed
const MAX_AGE: f64 = SOLAR_WIND_MIN_VELOCITY * 2.0 * 32.0 * (GRID_SIZE as f64 + 150.0) * 10.0;
const WIND_SPAWN_DISTANCE: f64 = 0.1 * (GRID_SIZE as f64);
const SPAWN_CHANCE: f64 = 0.5;

fn do_tick(tick: u32, indiv: &mut Individual) {
    // print!("Tick: {} ", tick);
    // println!("Particles: {:?}", particles.len());

    let mut rng = rand::rng();

    //Done
    if tick % (1 * SOLAR_WIND_TICK_MULTIPLIER) == 0 {
        tick_1_move_solar_wind(&mut indiv.particles)
    }

    if tick % (8 * SOLAR_WIND_TICK_MULTIPLIER) == 0 {
        tick_8_solar_wind_collisions(tick, &mut indiv.particles, &mut indiv.boxes);
    }

    //Done
    if tick % (SOLAR_WIND_DEFLECTION_TICK_INTERVAL * SOLAR_WIND_TICK_MULTIPLIER) == 0 {
        tick_solar_wind_deflection(&mut indiv.particles, &indiv.rods);
    }

    //Done - 9
    if tick % (9 * SOLAR_WIND_TICK_MULTIPLIER) == 0 {
        if rng.random_range(0.0..1.0) < SPAWN_CHANCE {
            for _ in 0..1 {
                spawn_solar_wind_particle(&mut indiv.particles);
                indiv.total_particles += 1;
            }
        }
    }
    if tick % 240 == 0 {
        tick_240_clean_up_cerys_solar_wind_particles(&mut indiv.particles)
    }
}

fn tick_1_move_solar_wind(particles: &mut Vec<Particle>) {
    for p in particles.iter_mut() {
        if !p.valid {
            continue;
        }
        p.x += p.vx;
        p.y += p.vy;
        p.age += 1;
    }
}

fn tick_8_solar_wind_collisions(tick: u32, particles: &mut Vec<Particle>, boxes: &mut Vec<Box>) {
    //this is the big one
    for p in particles.iter_mut() {
        let found_boxes: Vec<&mut Box> = boxes
            .iter_mut()
            .filter(|b| {
                let dx = p.x - b.x;
                let dy = p.y - b.y;
                let d2 = dx * dx + dy * dy;
                d2 < (0.75 * 0.75)
            })
            .collect();
        for b in found_boxes {
            if (!particle_is_in_cooldown(tick, p))
                || ((p.last_checked_inv > 0) && (p.last_checked_inv != b.id))
            {
                p.irradiation_tick = tick;
                p.last_checked_inv = b.id;
                b.count += 1;
                b.last_hit = tick;
            }
        }
    }
}

fn particle_is_in_cooldown(tick: u32, p: &mut Particle) -> bool {
    if p.irradiation_tick > 0 {
        return false;
    }

    let v2 = p.vx * p.vx + p.vy * p.vy;
    let speed = v2.sqrt();

    let cooldown_time_1 = (COOLDOWN_DISTANCE / speed).round() as u32;
    let cooldown_time_2 = COOLDOWN_TICKS as u32;

    if (tick > p.irradiation_tick + cooldown_time_1)
        || (tick > p.irradiation_tick + cooldown_time_2)
    {
        p.irradiation_tick = 0;
        p.last_checked_inv = 0;
        return false;
    }

    return true;
}

fn tick_solar_wind_deflection(particles: &mut Vec<Particle>, rods: &Vec<Rod>) {
    let mut rng = rand::rng();

    for p in particles.iter_mut() {
        for r in rods.iter() {
            let mut dx = p.x - r.x;
            let mut dy = p.y - r.y;
            let mut d2 = dx * dx + dy * dy;

            if d2 == 0.0 {
                let random_angle = rng.random::<f64>() * 2.0 * std::f64::consts::PI;
                dx = MIN_ELECTROMAGNETIC_DISTANCE * random_angle.cos();
                dy = MIN_ELECTROMAGNETIC_DISTANCE * random_angle.sin();
            } else if d2 < MIN_ELECTROMAGNETIC_DISTANCE * MIN_ELECTROMAGNETIC_DISTANCE {
                let scale = MIN_ELECTROMAGNETIC_DISTANCE / d2.sqrt();
                dx *= scale;
                dy *= scale;
                d2 = dx * dx + dy * dy;
            }

            if (d2 < ROD_MAX_RANGE_SQUARED) && (r.polarity_bool != 0.0) {
                let deflection = r.polarity_bool
                    * ROD_DEFLECTION_STRENGTH
                    * (SOLAR_WIND_DEFLECTION_TICK_INTERVAL as f64)
                    / 60.0;
                let dvx = dx / (d2.powf(7.0 / 4.0)) * deflection;
                let dvy = dy / (d2.powf(7.0 / 4.0)) * deflection;
                p.vx += dvx;
                p.vy += dvy;
            }
        }
    }
}

fn spawn_solar_wind_particle(particles: &mut Vec<Particle>) {
    let mut rng = rand::rng();

    let particle = Particle {
        //x: -(WIND_SPAWN_DISTANCE - rng.random_range(0.0..10.0)),
        x: -WIND_SPAWN_DISTANCE * rng.random::<f64>(),
        // y: rng.random_range(0.0..= (2.0 * (MOON_RADIUS + 8.0))),
        y: rng.random_range(-(GRID_SIZE as f64)-WIND_SPAWN_DISTANCE..=(GRID_SIZE as f64)+WIND_SPAWN_DISTANCE ),
        vx: SOLAR_WIND_MIN_VELOCITY + rng.random_range(0.0..1.0) * 0.05,
        vy: 0.3 * (rng.random::<f64>() - 0.5).powi(3),
        age: 0,
        valid: true,
        irradiation_tick: 0,
        last_checked_inv: 0,
    };
    particles.push(particle);
}

fn tick_240_clean_up_cerys_solar_wind_particles(particles: &mut Vec<Particle>) {
    for p in particles.iter_mut() {
        if p.age > MAX_AGE as u32 {
            p.valid = false;
        }
        // If they left the sim
        // else if (p.x.abs() > WIND_SPAWN_DISTANCE + 5.0) || (p.y.abs() > WIND_SPAWN_DISTANCE + 5.0)
        else if (p.x < -WIND_SPAWN_DISTANCE)
            || (p.x > (GRID_SIZE as f64) + WIND_SPAWN_DISTANCE)
            || (p.y < -WIND_SPAWN_DISTANCE)
            || (p.y > (GRID_SIZE as f64) + WIND_SPAWN_DISTANCE)
        {
            p.valid = false;
        }
    }

    particles.retain(|p| p.valid);
}

fn draw(window: &mut Window, tick: u32, indiv: &mut Individual) {
    let mut buffer = vec![0u32; WIDTH * HEIGHT];

    // Draw moon as a white circle, using grid_coord_is_on_moon
    // for x in 0..GRID_SIZE {
    //     for y in 0..GRID_SIZE {
    //         if grid_coord_is_on_moon(x, y) {
    //             let sx = grid_to_screenx(x);
    //             let sy = grid_to_screenx(y);
    //             for dx in 0..SCALE {
    //                 for dy in 0..SCALE {
    //                     if sx + dx < WIDTH && sy + dy < HEIGHT {
    //                         buffer[(sy + dy) * WIDTH + (sx + dx)] = 0x606060;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // Draw grid as gray lines
    for x in 0..WIDTH {
        for y in 0..HEIGHT {
            if x % SCALE == 0 || y % SCALE == 0 {
                buffer[y * WIDTH + x] = 0x303030;
            }
        }
    }

    // Draw rods as red squares
    for rod in indiv.rods.iter() {
        let color = if rod.polarity_bool > 0.0 {
            0x0000FF
        } else {
            0xFF0000
        };
        let rx = rod.gridx * SCALE;
        let ry = rod.gridy * SCALE;
        let margin = (SCALE as f64 * 0.1) as usize;
        if rx < WIDTH && ry < HEIGHT {
            for sx in margin..(SCALE * 2 - margin) {
                for sy in margin..(SCALE * 2 - margin) {
                    if rx + sx < WIDTH && ry + sy < HEIGHT {
                        buffer[(ry + sy) * WIDTH + (rx + sx)] = color;
                    }
                }
            }
        }
    }

    // Draw boxes as green squares 1x1 with margin
    for b in indiv.boxes.iter() {
        let color = {
            if b.last_hit > 0 {
                let elapsed_ticks = tick.saturating_sub(b.last_hit);
                let lerp_factor = 1.0 - (elapsed_ticks as f64 / 60.0).min(1.0);
                let r = 0x00;
                let g = (0x88 + ((0xFF - 0x88) as f64 * lerp_factor) as u32) & 0xFF;
                let b = 0x00;
                (r << 16) | (g << 8) | b
            } else {
                0x008800
            }
        };
        let rx = b.gridx * SCALE;
        let ry = b.gridy * SCALE;
        let margin = (SCALE as f64 * 0.1) as usize;
        if rx < WIDTH && ry < HEIGHT {
            for sx in margin..(SCALE - margin) {
                for sy in margin..(SCALE - margin) {
                    if rx + sx < WIDTH && ry + sy < HEIGHT {
                        buffer[(ry + sy) * WIDTH + (rx + sx)] = color;
                    }
                }
            }
        }
    }

    // Draw particles as white dots
    const PINK: u32 = 0xFF00FF; // 50% transparent pink
    for p in indiv.particles.iter() {
        let px = (p.x * SCALE as f64) as isize;
        let py = (p.y * SCALE as f64) as isize;
        if px >= 0 && py >= 0 && (px as usize) < WIDTH && (py as usize) < HEIGHT {
            for dx in -(PARTICLE_RADIUS as isize)..=(PARTICLE_RADIUS as isize) {
                for dy in -(PARTICLE_RADIUS as isize)..=(PARTICLE_RADIUS as isize) {
                    if dx * dx + dy * dy <= (PARTICLE_RADIUS * PARTICLE_RADIUS) as isize {
                        let nx = px + dx;
                        let ny = py + dy;
                        if nx >= 0 && ny >= 0 && (nx as usize) < WIDTH && (ny as usize) < HEIGHT {
                            buffer[(ny as usize) * WIDTH + (nx as usize)] = PINK;
                        }
                    }
                }
            }
        }
    }

    window.update_with_buffer(&buffer, WIDTH, HEIGHT).unwrap();
}

#[derive(Debug, Clone)]
struct Particle {
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    age: u32,
    valid: bool,
    irradiation_tick: u32,
    last_checked_inv: u32,
}

#[derive(Debug, Clone)]
struct Rod {
    gridx: usize,
    gridy: usize,
    x: f64,
    y: f64,
    polarity: f64,
    polarity_bool: f64,
}

impl Rod {
    fn new(gridx: usize, gridy: usize, polarity: f64) -> Self {
        Rod {
            gridx: gridx,
            gridy: gridy,
            x: (gridx + 1) as f64,
            y: (gridy + 1) as f64,
            polarity: polarity,
            polarity_bool: if polarity > 0.0 { 1.0 } else { -1.0 },
        }
    }
}

#[derive(Debug, Clone)]
struct Box {
    gridx: usize,
    gridy: usize,
    x: f64,
    y: f64,
    count: u32,
    id: u32,
    last_hit: u32,
}

impl Box {
    fn new(id: u32, gridx: usize, gridy: usize) -> Self {
        Box {
            gridx: gridx,
            gridy: gridy,
            x: (gridx as f64) + 0.5,
            y: (gridy as f64) + 0.5,
            count: 0,
            id: id,
            last_hit: 0,
        }
    }
}

#[derive(Debug, Clone)]
struct Individual {
    id: u32,
    max_rods: usize,
    max_boxes: usize,
    rods: Vec<Rod>,
    boxes: Vec<Box>,
    sim_size: usize,
    particles: Vec<Particle>,
    total_particles: u32,
}

impl Individual {
    fn grid_location_taken(rods: &Vec<Rod>, boxes: &Vec<Box>, x: usize, y: usize) -> bool {
        let taken_by_boxes = boxes.iter().any(|b| b.gridx == x && b.gridy == y);
        let taken_by_rods = rods
            .iter()
            .any(|r| (r.gridx == x || r.gridx + 1 == x) && (r.gridy == y || r.gridy + 1 == y));
        taken_by_boxes || taken_by_rods
    }

    fn random_placement(&mut self) {
        let mut rng = rand::rng();

        while self.boxes.len() < self.max_boxes {
            let x = rng.random_range(0..self.sim_size - 1);
            let y = rng.random_range(0..self.sim_size - 1);

            if !Self::grid_location_taken(&self.rods, &self.boxes, x, y) {
                self.boxes.push(Box::new(rng.random::<u32>(), x, y));
            }
        }

        while self.rods.len() < self.max_rods {
            let x = rng.random_range(1..self.sim_size - 1);
            let y = rng.random_range(1..self.sim_size - 1);

            if !(Self::grid_location_taken(&self.rods, &self.boxes, x, y)
                || Self::grid_location_taken(&self.rods, &self.boxes, x + 1, y)
                || Self::grid_location_taken(&self.rods, &self.boxes, x, y + 1)
                || Self::grid_location_taken(&self.rods, &self.boxes, x + 1, y + 1))
            {
                self.rods.push(Rod::new(x, y, rng.random::<f64>() - 0.5));
            }
        }
    }

    fn new(sim_size: usize, max_rods: usize, max_boxes: usize) -> Self {
        let mut rng = rand::rng();
        let mut indiv = Individual {
            id: rng.random::<u32>(),
            max_rods,
            max_boxes,
            rods: Vec::new(),
            boxes: Vec::new(),
            sim_size,
            particles: Vec::new(),
            total_particles: 0,
        };

        indiv.random_placement();

        indiv
    }

    fn serialize(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("{}\n", self.sim_size));
        s.push_str(&format!("{}\n", self.rods.len()));
        for r in self.rods.iter() {
            s.push_str(&format!("{} {} {} {}\n", r.gridx, r.gridy, r.polarity, r.polarity_bool));
        }
        s.push_str(&format!("{}\n", self.boxes.len()));
        for b in self.boxes.iter() {
            s.push_str(&format!("{} {}\n", b.gridx, b.gridy));
        }
        s
    }

    fn from_string(s: &str) -> Self {
        let mut lines = s.lines();
        let sim_size = lines.next().unwrap().parse::<usize>().unwrap();
        let rod_count = lines.next().unwrap().parse::<usize>().unwrap();
        let mut rods = Vec::new();
        for _ in 0..rod_count {
            let mut parts = lines.next().unwrap().split_whitespace();
            let gridx = parts.next().unwrap().parse::<usize>().unwrap();
            let gridy = parts.next().unwrap().parse::<usize>().unwrap();
            let polarity = parts.next().unwrap().parse::<f64>().unwrap();
            let polarity_bool = parts.next().unwrap().parse::<f64>().unwrap();
            rods.push(Rod {
                gridx,
                gridy,
                x: (gridx + 1) as f64,
                y: (gridy + 1) as f64,
                polarity,
                polarity_bool,
            });
        }
        let box_count = lines.next().unwrap().parse::<usize>().unwrap();
        let mut boxes = Vec::new();
        for _ in 0..box_count {
            let mut parts = lines.next().unwrap().split_whitespace();
            let gridx = parts.next().unwrap().parse::<usize>().unwrap();
            let gridy = parts.next().unwrap().parse::<usize>().unwrap();
            boxes.push(Box::new(0, gridx, gridy));
        }
        let mut rng = rand::rng();
        Individual {
            id: rng.random::<u32>(),
            max_rods: rods.len(),
            max_boxes: boxes.len(),
            rods,
            boxes,
            sim_size,
            particles: Vec::new(),
            total_particles: 0,
        }
    }

    fn hits_per_particle_per_box(&self) -> f64 {
        self.boxes.iter().map(|b| b.count).sum::<u32>() as f64 / (self.boxes.len() as f64 * self.total_particles as f64)
    }


    fn fitness(&self) -> f64 {
        let boxes_adjustment =  1.0 - (self.boxes.len() as f64 / MAX_BOXES as f64);
        let rods_adjustment = 1.0 - (self.rods.len() as f64 / MAX_RODS as f64);
        let effiency = rods_adjustment * boxes_adjustment;
        // println!("Boxes {} Rods {} Effiency: {}", self.boxes.len(), self.rods.len(), effiency);

        self.hits_per_particle_per_box()
    }

    fn reset_sim(&mut self) {
        self.particles.clear();
        self.total_particles = 0;
        self.boxes.iter_mut().for_each(|b| {
            b.count = 0;
            b.last_hit = 0
        });
    }

    fn random_mutate(&mut self, alpha: f64) {
        let mut rng = rand::rng();

        if alpha < rng.random::<f64>(){
            self.max_rods = (self.max_rods as f64 + (4.0 * (rng.random::<f64>() - 0.5)))
                .round() as usize;
            self.max_rods = self.max_rods.clamp(1, MAX_RODS);
        }
        if alpha < rng.random::<f64>(){
            self.max_boxes = (self.max_boxes as f64 + (4.0 * (rng.random::<f64>() - 0.5)))
                .round() as usize;
            self.max_boxes = self.max_boxes.clamp(1, MAX_BOXES);
        }

        if self.boxes.len() > self.max_boxes {
            while self.boxes.len() > self.max_boxes {
                let index = rng.random_range(0..self.boxes.len());
                self.boxes.remove(index);
            }
        }

        if self.rods.len() > self.max_rods {
            while self.rods.len() > self.max_rods {
                let index = rng.random_range(0..self.rods.len());
                self.rods.remove(index);
            }
        }

        self.random_placement();

        self.rods.iter_mut().for_each(|r| {
            r.polarity += alpha * 2.0 * (rng.random::<f64>() - 0.5);
            r.polarity_bool = if r.polarity > 0.0 { 1.0 } else { -1.0 };
        });
    }

    fn breed_individuals(i1: Individual, i2: Individual) -> Individual {
        
        let mut rng = rand::rng();
    
        let mut max_rods = (i1.max_rods + i2.max_rods / 2).clamp(1, MAX_RODS);
        let mut max_boxes = (i1.max_boxes + i2.max_boxes / 2).clamp(1, MAX_BOXES);
        
        let mut boxes: Vec<Box> = Vec::new();
        let mut rods: Vec<Rod> = Vec::new();

        let mut tries = 0;
        while (boxes.len() < max_boxes) && (tries < 1000) {
            tries += 1;
            let b = if rng.random::<f64>() < 0.5 {
                i1.boxes[rng.random_range(0..i1.boxes.len())].clone()
            } else {
                i2.boxes[rng.random_range(0..i2.boxes.len())].clone()
            };
            if !Self::grid_location_taken(&rods, &boxes, b.gridx, b.gridy) {
                boxes.push(b);
            }
        }
        max_boxes = boxes.len();
        
    
        tries = 0;
        while (rods.len() < max_rods) && (tries < 1000) {
            tries += 1;
            let rod = if rng.random::<f64>() < 0.5 {
                i1.rods[rng.random_range(0..i1.rods.len())].clone()
            } else {
                i2.rods[rng.random_range(0..i2.rods.len())].clone()
            };
            if !(Self::grid_location_taken(&rods, &boxes, rod.gridx, rod.gridy)
                || Self::grid_location_taken(&rods, &boxes, rod.gridx + 1, rod.gridy)
                || Self::grid_location_taken(&rods, &boxes, rod.gridx, rod.gridy + 1)
                || Self::grid_location_taken(&rods, &boxes, rod.gridx + 1, rod.gridy + 1))
            {
                rods.push(rod);
            }
        }
        max_rods = rods.len();
        
        Individual {
            id: rng.random::<u32>(),
            max_rods: max_rods.clamp(1, MAX_RODS),
            max_boxes: max_boxes.clamp(1, MAX_BOXES),
            rods: rods,
            boxes: boxes,
            sim_size: i1.sim_size,
            particles: Vec::new(),
            total_particles: 0,
        }
    
    }
}



fn main() {
    let mut rng = rand::rng();

    let mut population: Vec<Individual> = (0..POPULATION_SIZE)
        .into_par_iter()
        .map(|_| Individual::new(GRID_SIZE, NUM_RODS, NUM_BOXES))
        .collect::<Vec<Individual>>();

    for gen_num in 1..=NUM_GENERATIONS {
        population
            .par_iter_mut()
            .progress_count(POPULATION_SIZE)
            .for_each(|indiv| {
                let mut tick = 0;

                while tick < 20_000 {
                    do_tick(tick, indiv);
                    tick += 1;
                }
            });

        //sort
        population.sort_by(|a, b| b.fitness().partial_cmp(&a.fitness()).unwrap());

        let average_fitness: f64 =
            population.iter().map(|indiv| indiv.fitness()).sum::<f64>() / population.len() as f64;
        //best
        let mut fittest: Individual = population[0].clone();
        
        //
        let percentiles = [0.01, 0.20, 0.50];

        let fitness_values: Vec<f64> = population.iter().map(|indiv| indiv.fitness()).collect();
        let hits_values: Vec<f64>= population.iter().map(|indiv| indiv.hits_per_particle_per_box()).collect();
        let box_counts: Vec<usize> = population.iter().map(|indiv| indiv.boxes.len()).collect();
        let rod_counts: Vec<usize> = population.iter().map(|indiv| indiv.rods.len()).collect();

        let get_percentile_value = |values: &Vec<f64>, percentile: f64| -> f64 {
            let index = (percentile * values.len() as f64).ceil() as usize - 1;
            values[index]
        };

        let get_percentile_count = |counts: &Vec<usize>, percentile: f64| -> usize {
            let index = (percentile * counts.len() as f64).ceil() as usize - 1;
            counts[index]
        };
        println!("Generation: {}", gen_num);
        println!("{:<15} {:<35} {:<15} {:<15} {:<15}", "Percentile", "Hits per Particle per Box", "Fitness", format!("Boxes\\{}", MAX_BOXES), format!("Rods\\{}", MAX_RODS));
        for &percentile in &percentiles {
            println!(
            "{:<15.0} {:<35.5} {:<15.5} {:<15} {:<15}",
            percentile * 100.0,
            get_percentile_value(&fitness_values, percentile),
            get_percentile_value(&hits_values, percentile),
            get_percentile_count(&box_counts, percentile),
            get_percentile_count(&rod_counts, percentile)
            );
        }
        println!("");

        

        if (gen_num % 10) == 0{

            // Save fittest individual to file
            let filename = format!("fit/fittest_gen_{}.txt", gen_num);
            let mut file = File::create(filename).unwrap();
            file.write_all(fittest.serialize().as_bytes()).unwrap();

            
            let mut window = Window::new(
                "Particle Simulation",
                WIDTH,
                HEIGHT,
                WindowOptions::default(),
            )
            .unwrap();
            window.set_target_fps(60);

            let mut tick: u32 = 0;
            while window.is_open() {
                do_tick(tick, &mut fittest);
                draw(&mut window, tick, &mut fittest);

                thread::sleep(FRAME_TIME);
                tick += 1;
            }
        }

        let ranks: Vec<_> = population
            .iter()
            .map(|indiv| {
                let rank = population.iter().position(|x| x.id == indiv.id).unwrap();
                (indiv.id, rank)
            })
            .collect();

        let surviving_population: Vec<_> = population
            .iter()
            .filter(|indiv| {
                let rank = ranks.iter().find(|&&(id, _)| id == indiv.id).unwrap().1;
                let survival_chance = 1.0 - (rank as f64 / POPULATION_SIZE as f64);
                rng.random::<f64>() <= survival_chance
            })
            .cloned()
            .collect();

        for indiv in population.iter_mut() {
            let rank = ranks.iter().find(|&&(id, _)| id == indiv.id).unwrap().1;
            let survival_chance = 1.0 - (rank as f64 / POPULATION_SIZE as f64);

            
            if rng.random::<f64>() > 0.2 {
                //Asexual reproduction
                if rng.random::<f64>() > survival_chance {
                    let random_survivor = surviving_population
                        [rng.random_range(0..surviving_population.len())]
                    .clone();
                    let mut new_indiv = random_survivor.clone();
                    if rng.random::<f64>() > 0.5 {
                        new_indiv.random_mutate(0.1);
                    }
                    *indiv = new_indiv;
                } else {
                    //Sexual reproduction
                    let parent1 = surviving_population
                        [rng.random_range(0..surviving_population.len())]
                    .clone();
                    let parent2 = surviving_population
                        [rng.random_range(0..surviving_population.len())]
                    .clone();
                    let mut child = Individual::breed_individuals(parent1, parent2);
                    child.random_mutate(0.02);
                    *indiv = child;
                }
            }
        }
    }
}
