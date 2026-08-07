# MHD Nanofluid Couette Flow: Interactive Simulator

An interactive web application for the steady-state analysis of magnetohydrodynamic (MHD)
nanofluid Couette flow with partial slip and convective cooling. It solves the coupled
momentum and energy equations live as you change parameters, visualises velocity,
temperature, entropy generation and the Bejan number, and includes a trained machine-learning
surrogate for instant prediction of the engineering outputs.

## About

The model describes a single-phase nanofluid between two parallel plates. The lower plate is
stationary and isothermal at ambient temperature; the upper plate moves with Navier partial
slip and exchanges heat convectively with the surroundings, under a uniform transverse
magnetic field. The app lets you explore how the flow, heat transfer and thermodynamic
irreversibility respond to the governing dimensionless groups.

### The science in brief

Dimensionless governing equations (steady state):

```
Momentum:  A1 * W'' - A2 * Ha^2 * W + G = 0
Energy:    A3 * theta'' + A1 * Pr * Ec * (W')^2 + A2 * Pr * Ec * Ha^2 * W^2 = 0
```

Boundary conditions:

```
eta = 0 :  W = 0,                 theta = 0           (stationary, isothermal lower plate)
eta = 1 :  W = Re - lambda*W'(1), theta'(1) = -Bi*theta(1)   (partial slip + convective cooling)
```

Engineering outputs:

```
Skin friction:  Cf = A1 * (dW/deta) at the wall
Nusselt number: Nu = -A3 * (dtheta/deta) at the wall     (SIGNED: sign gives heat-flux direction)
Entropy gen.:   Ns = A3*(theta')^2 + A1*Ec*Pr*(W')^2 + A2*Pr*Ec*Ha^2*W^2
Bejan number:   Be = N1 / (N1 + N2 + N3)
```

Nu is reported as a signed quantity: the magnitude is the heat-transfer rate and the sign is
the direction of wall heat flux. A notable result is the lower-wall Nu(0) sign reversal near
Re around 3 to 4, where the stationary wall switches from heat sink to heat source.

Nanofluid effective properties use standard models: Brinkman (viscosity), Maxwell (electrical
conductivity) and Maxwell-Garnett (thermal conductivity). The ratios A1, A2, A3 in the
equations are derived from these correlations.

### Dimensionless groups

| Symbol | Name | Meaning |
|--------|------|---------|
| Ha | Hartmann number | Magnetic to viscous force ratio |
| Re | Reynolds number | Upper-plate motion |
| Pr | Prandtl number | Momentum to thermal diffusivity |
| Ec | Eckert number | Viscous dissipation |
| Bi | Biot number | Convective to conductive resistance |
| lambda | Slip parameter | Navier partial slip at the moving wall |
| G | Pressure gradient | Axial driving pressure |
| A1, A2, A3 | Property ratios | Viscosity, electrical and thermal conductivity of nanofluid vs base fluid |

## Numerical method (honest note)

The published research results are produced with the Spectral Quasi-linearisation Method
(SQLM), which reaches near machine precision. This web application uses a second-order Finite
Difference Method (FDM) instead, chosen for fast, interactive, in-browser solving. The two
methods are validated against the same four exact analytical limiting cases; the FDM errors
are larger than SQLM's, as expected, but small and convergent. The Validation tab computes
these comparisons live.

## Machine-learning surrogate

The AI Lab tab runs a trained scikit-learn MLPRegressor (architecture 11 to 256 to 256 to 128
to 64 to 9, ReLU) that learned the SQLM solver's input to output map. The exported weights and
scalers are loaded at runtime from `public/mlp_model.json` and the forward pass runs entirely
in the browser. Mean test accuracy is R^2 = 0.9958 against the SQLM ground truth. It predicts
all nine engineering outputs; the displayed value is the prediction next to the live FDM solve.

## Tech stack

- React (Create React App)
- Recharts for plotting
- lucide-react for icons
- Supabase JavaScript client (optional, for an opt-in community dataset feature)
- Model trained offline with scikit-learn; only the exported weights ship in the app

## Getting started

Prerequisites: Node.js LTS (18 or newer) and npm.

```bash
# install dependencies
npm install

# run the development server (http://localhost:3000)
npm start

# create an optimised production build in ./build
npm run build
```

## Environment variables (optional)

The app runs fully offline without any configuration. Supabase is only used for an optional,
opt-in community dataset feature, which is disabled by default. To enable it, provide your own
project values:

```
REACT_APP_SUPABASE_URL=your-project-url
REACT_APP_SUPABASE_ANON_KEY=your-anon-key
```

Only ever use a Supabase anon or publishable key in client code, never a service-role key, and
enable Row-Level Security on any tables the app touches.

## Project structure

```
public/
  index.html
  mlp_model.json      # exported ML surrogate weights and scalers
  images/             # key research figures
  videos/             # explainer videos
src/
  App.js              # solver, UI, tabs, ML forward pass
  supabaseClient.js   # optional Supabase client
  hooks/              # data hook for the optional community feature
  index.js, index.css
```

## Deployment

The production build in `./build` is a static site and can be hosted on any static host.

## Authors and acknowledgements

- Candidate: Mr. S.I. Mosala
- Supervisors: Prof. O.D. Makinde and Prof. A. Muronga
- Institution: Nelson Mandela University

## License

For educational and research purposes.
