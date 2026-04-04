***

# LQR and Predictive Control of the "Ball and Beam" Model ⚙️🎱

This repository contains the MATLAB code and documentation for the modeling and optimal control of a **Ball and Beam** system. This work was developed for the *Optimal Control* course in the Robotics and Automation Engineering degree program at the University of Calabria. 
Author: Giuseppe Coppola.

## 📖 Project Description

The Ball and Beam system is a classic example that represents an unstable and nonlinear system, ideal for putting advanced regulation techniques into practice. The goal of the control system is to regulate the angle of the beam (via the torque of an actuator) in order to maintain the ball in a desired position (typically near the center of the beam), compensating for physical perturbations.

The project explores various control strategies, starting from fundamental linear control up to robust Model Predictive Control (MPC) techniques based on Linear Matrix Inequalities (LMI) to handle model uncertainties and strict actuation constraints.

## 🚀 Implemented Control Strategies

By running the code, the following conditions and strategies will be simulated sequentially:

1. **Free Response:** Analysis of the open-loop system (the system is naturally unstable).
2. **LQR Control (Infinite Horizon):** Zero regulation using Dynamic Programming to minimize a quadratic performance index.
3. **LQR Control (Finite Horizon):** Application of a time-varying sequence of feedback gains.
4. **Robust Unconstrained RHC:** Unconstrained robust Receding Horizon Control, based on a polytopic description of the linearized model's uncertainties.
5. **Robust Constrained Input RHC:** RHC with strict saturation constraints imposed on the motor input ($u_{max}$).
6. **Robust Constrained Input and Output RHC:** Extension of the RHC controller to impose *component-wise* constraints not only on the input but directly on the state variables (position, velocity, angle).

## 🛠️ System Requirements

To properly run all scripts, you need to have the following installed:
* **MATLAB** (Tested on recent environments).
* **Control System Toolbox**.
* **Symbolic Math Toolbox** (used to derive Jacobians for automatic linearization).
* **Optimization Toolbox** (used in `build_pol.m` via `linprog`).
* **YALMIP**: Advanced parser for optimization in MATLAB.
* **SeDuMi**: Semidefinite Programming (SDP) solver configured in YALMIP to solve the LMI problems present in the `MPC` class.

## 📂 Repository Structure

* `Main.m`: The main script. Initializes system parameters (mass, moment of inertia, etc.), defines the dynamical model using Lagrangian Mechanics, and runs all simulations progressively.
* `OptControl.m`: MATLAB class dedicated to linearizing the model (first-order Taylor expansion) and synthesizing the infinite and finite horizon LQR controllers.
* `MPC.m`: MATLAB class that implements the robust RHC controllers. It sets up the optimization problems and solves them as LMIs using YALMIP.
* `build_pol.m`: Function dedicated to the automatic construction of the uncertain model. It samples the state space and builds the exact convex-hull for the polytopic representation, avoiding redundant vertices.
* `ell.m`: Auxiliary class that handles instantiating, validating, and plotting the ellipsoids relative to the Robust Positively Invariant (RPI) sets, which are used to verify asymptotic stability.
* `Ball_Beam_OC.pdf`: Complete technical report. Contains the full mathematical derivation of the Lagrangian, analysis of the graphical results, and the theoretical formulation behind the code.

## 💻 How to Use

1. Ensure **YALMIP** and **SeDuMi** are correctly installed and added to your MATLAB environment's *Path*.
2. Clone this repository to a local folder.
3. Open MATLAB and set your working directory to the folder containing the files.
4. Run the `Main.m` script.
5. The script is interactive: it will pause after each simulation. Click inside the MATLAB *Command Window* and press `Enter` to move to the next control algorithm. All state evolution, input, and invariant set (RPI) plots will be generated automatically and sequentially.
