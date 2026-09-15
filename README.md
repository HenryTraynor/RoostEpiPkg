# RoostEpiPkg

## Spatial movement and epidemic modeling in structured wildlife populations

`RoostEpiPkg` is an R package for simulating the interaction between spatial movement and epidemic dynamics in structured wildlife populations.

The model combines:

- discrete populations located at spatially distributed roosts
- SIR or SIRS epidemiological dynamics within roosts
- continuous spatial movement across a landscape
- distance-dependent movement
- spatially varying diffusion
- population- and capacity-dependent roost attraction
- capacity-constrained assignment of individuals back to roosts

The package was developed to investigate how spatial structure and movement influence epidemic dynamics in wildlife populations.

---

## Model overview

Each simulation day consists of three major processes:

```text
                    DAY
                     │
                     ▼
             ┌───────────────┐
             │  Epidemiology │
             │   in roosts   │
             └───────┬───────┘
                     │
                     ▼
             ┌───────────────┐
             │   Nighttime   │
             │    movement   │
             │  continuous   │
             │    diffusion  │
             └───────┬───────┘
                     │
                     ▼
             ┌───────────────┐
             │ Dawn return & │
             │ roost assign. │
             └───────┬───────┘
                     │
                     ▼
                  NEXT DAY
```

### Daytime epidemiology

Individuals are assigned to discrete roosts during the day. Disease transmission occurs within each roost using either an SIR or SIRS model.

For the SIR model:

```text
S → I → R
```

For the SIRS model:

```text
S → I → R → S
```

Transmission and recovery are implemented as stochastic transitions.

### Nighttime movement

At night, individuals leave the discrete roost representation and are represented as continuous spatial densities across a two-dimensional landscape.

Movement is modeled as diffusion with spatially varying diffusion coefficients.

The movement kernel uses exponential distance decay:

$$
K(d) = e^{-\lambda d}
$$

where:

- $d$ is distance
- $\lambda$ controls the strength of distance decay

Larger values of $\lambda$ produce stronger distance decay.

### Roost connectivity

Connectivity between roosts is calculated using the same distance-decay function:

$$
C_{jk} = e^{-\lambda d_{jk}}
$$

where $d_{jk}$ is the distance between roosts $j$ and $k$.

Roost isolation is derived from total connectivity and influences the diffusion coefficient assigned to the roost and surrounding landscape.

### Dawn return

At dawn, continuous spatial densities are converted back into discrete roost populations.

The probability that an individual in cell $i$ returns to roost $j$ is based on:

$$
W_{ij} = K(d_{ij}) A_j
$$

where $A_j$ represents population- and capacity-dependent attraction to roost $j$.

Roost capacity is explicitly enforced during assignment.

---

## Installation

The development version can be installed from GitHub using:

```r
install.packages("remotes")

remotes::install_github("YOUR-GITHUB-USERNAME/RoostEpiPkg")
```

Once the package is released, the installation instructions will be updated accordingly.

---

## Quick start

Load the package:

```r
library(RoostEpiPkg)
```

### 1. Create roosts

```r
roosts <- create_roosts(
  num_roosts = 25,
  num_clusters = 5,
  sd = 0.05,
  avg_N_max = 70,
  min_N_max = 10
)
```

### 2. Initialize populations

```r
roosts <- initialize_roost_populations(
  roosts,
  initial_occupancy = 0.25,
  initial_infected = 0.25
)
```

### 3. Create the spatial grid

```r
spatial <- create_spatial_params(
  x_min = 0,
  x_max = 1,
  y_min = 0,
  y_max = 1,
  nx = 100,
  ny = 100
)
```

### 4. Set up movement

```r
movement <- setup_movement(
  roost_info = roosts,
  spatial_params = spatial,
  lambda = 5,
  sigma = 0.01,
  D_base = 0.0005,
  gamma_D = 1,
  alpha = 1,
  epsilon = 0.01
)
```

### 5. Set up the epidemic model

```r
epidemic <- create_epidemic_params(
  beta = 0.5,
  gamma = 0.1,
  model = "SIR"
)
```

### 6. Run the simulation

```r
sim <- run_simulation(
  roost_info = roosts,
  spatial_params = spatial,
  movement_params = movement,
  epidemic_params = epidemic,
  n_days = 30
)
```

The resulting object contains the model state for every simulated day.

---

## Model parameters

### Movement

| Parameter | Description |
|---|---|
| `lambda` | Strength of exponential distance decay |
| `sigma` | Spatial spread of the initial roost density |
| `D_base` | Baseline diffusion coefficient |
| `gamma_D` | Effect of roost isolation on diffusion |
| `alpha` | Strength of population-based roost attraction |
| `epsilon` | Small value preventing zero attraction |
| `times_night` | Time sequence used for nighttime diffusion |

### Epidemiology

| Parameter | Description |
|---|---|
| `beta` | Transmission rate |
| `gamma` | Recovery rate |
| `omega` | Loss-of-immunity rate in the SIRS model |
| `model` | `"SIR"` or `"SIRS"` |

### Roosts

| Parameter | Description |
|---|---|
| `num_roosts` | Number of roosts |
| `num_clusters` | Number of spatial roost clusters |
| `sd` | Spatial spread around cluster centers |
| `avg_N_max` | Mean roost carrying capacity |
| `min_N_max` | Minimum roost carrying capacity |

---

## Model outputs

`run_simulation()` returns a list containing the state of the system for each simulated day.

For example:

```r
sim$day_0
sim$day_1
sim$day_2
```

Each day contains the roost populations:

```r
sim$day_1$roosts
```

and, for simulated movement, the nighttime spatial densities:

```r
sim$day_1$night$u_S
sim$day_1$night$u_I
sim$day_1$night$u_R
```

The dawn return probabilities are also retained:

```r
sim$day_1$roost_probability
```

---

## Population conservation

The model explicitly maintains several population-level constraints.

For every roost:

$$
N_j = S_j + I_j + R_j
$$

and:

$$
N_j \leq N_{max,j}
$$

The total population is conserved during nighttime diffusion and during the dawn reassignment process.

These properties are tested automatically using `testthat`.

---

## Documentation

A detailed description of the model is available in the package vignette:

```r
vignette(
  "continuous_model",
  package = "RoostEpiPkg"
)
```

The vignette describes the model structure, movement process, epidemiological processes, dawn return mechanism, and simulation workflow.

---

## Testing

The package contains automated tests for:

- roost creation and initialization
- spatial grid construction
- SIR and SIRS dynamics
- daytime epidemiology
- movement and diffusion
- population conservation
- dawn roost assignment
- roost capacity constraints
- complete multi-day simulations

Run the test suite with:

```r
devtools::test()
```

The package is also checked using:

```r
devtools::check()
```

---

## Development status

`RoostEpiPkg` is currently under active development.

The model structure and implementation are being developed for research applications involving spatial movement and infectious disease dynamics in wildlife populations.

---

## License

MIT
