# MHD Nanofluid Couette Flow - Interactive Simulation

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![React](https://img.shields.io/badge/React-18.2.0-61DAFB.svg)
![Status](https://img.shields.io/badge/status-Research%20Project-green.svg)

## 🌊 About This Project

This is an interactive web application for visualizing and simulating **Magnetohydrodynamic (MHD) Nanofluid Couette Flow**. The application solves the coupled momentum and energy equations in real-time, allowing users to explore how different parameters affect velocity, temperature, and entropy generation.

### Research Information
- **Project:** Thermal and Magnetohydrodynamic Analysis of Nanofluid Couette Flow
- **Candidate:** Mr. S.I. Mosala
- **Supervisor:** Prof. O.D. Makinde
- **Institution:** Nelson Mandela University
- **Date:** December 2025

---

## 📋 Table of Contents

1. [Prerequisites](#prerequisites)
2. [Step-by-Step Setup Guide](#step-by-step-setup-guide)
3. [Running the Application](#running-the-application)
4. [Deploying to Netlify via GitHub](#deploying-to-netlify-via-github)
5. [Adding Your Own Content](#adding-your-own-content)
6. [Project Structure](#project-structure)
7. [Physics Equations](#physics-equations)
8. [Troubleshooting](#troubleshooting)

---

## 🔧 Prerequisites

Before starting, you need to install these tools on your computer:

### 1. Install Node.js

Node.js is required to run the application.

1. Go to: https://nodejs.org/
2. Download the **LTS (Long Term Support)** version
3. Run the installer and follow the prompts
4. To verify installation, open Terminal (Mac) or Command Prompt (Windows) and type:
   ```bash
   node --version
   npm --version
   ```
   You should see version numbers (e.g., v18.x.x and 9.x.x)

### 2. Install Visual Studio Code

1. Go to: https://code.visualstudio.com/
2. Download and install for your operating system
3. Open VS Code after installation

### 3. Install Git

Git is needed to upload your code to GitHub.

1. Go to: https://git-scm.com/downloads
2. Download and install for your operating system
3. Verify installation:
   ```bash
   git --version
   ```

### 4. Create Accounts

- **GitHub Account:** https://github.com (free)
- **Netlify Account:** https://netlify.com (free, can sign up with GitHub)

---

## 📝 Step-by-Step Setup Guide

### Step 1: Create Your Project Folder

1. Create a new folder on your computer called `mhd-nanofluid-app`
2. You can put it anywhere (e.g., Desktop or Documents)

### Step 2: Copy the Project Files

Copy all the files I provided into your `mhd-nanofluid-app` folder. Your folder structure should look like this:

```
mhd-nanofluid-app/
├── public/
│   └── index.html
├── src/
│   ├── App.js
│   ├── index.js
│   └── index.css
├── package.json
└── README.md
```

### Step 3: Open in VS Code

1. Open VS Code
2. Click **File** → **Open Folder**
3. Navigate to and select your `mhd-nanofluid-app` folder
4. Click **Open**

### Step 4: Open the Terminal in VS Code

1. In VS Code, click **Terminal** in the top menu
2. Click **New Terminal**
3. A terminal panel will open at the bottom of VS Code

### Step 5: Install Dependencies

In the VS Code terminal, type the following command and press Enter:

```bash
npm install
```

⏳ Wait for the installation to complete (this may take 1-3 minutes). You'll see a progress bar and eventually a success message.

**Note:** You might see some warnings - this is normal and can be ignored as long as there are no errors.

---

## 🚀 Running the Application

### Start the Development Server

In the VS Code terminal, type:

```bash
npm start
```

### What Happens Next

1. The terminal will show "Compiled successfully!"
2. Your web browser will automatically open
3. The app will be running at: `http://localhost:3000`

### Stopping the Server

To stop the application, press `Ctrl + C` (Windows/Linux) or `Cmd + C` (Mac) in the terminal.

---

## 🌐 Deploying to Netlify via GitHub

### Part A: Upload to GitHub

#### Step 1: Create a GitHub Repository

1. Go to https://github.com
2. Click the **+** button (top right) → **New repository**
3. Name it: `mhd-nanofluid-simulation`
4. Set to **Public**
5. Do NOT initialize with README (uncheck this option)
6. Click **Create repository**

#### Step 2: Initialize Git in Your Project

In VS Code terminal, run these commands one at a time:

```bash
# Initialize git
git init

# Add all files
git add .

# Create first commit
git commit -m "Initial commit: MHD Nanofluid Couette Flow Simulation"
```

#### Step 3: Connect to GitHub

Replace `YOUR-USERNAME` with your actual GitHub username:

```bash
# Add your GitHub repository as remote
git remote add origin https://github.com/YOUR-USERNAME/mhd-nanofluid-simulation.git

# Rename branch to main
git branch -M main

# Push to GitHub
git push -u origin main
```

**Note:** If this is your first time using Git, you may be prompted to log in to GitHub.

### Part B: Deploy to Netlify

#### Step 1: Connect Netlify to GitHub

1. Go to https://netlify.com
2. Click **Sign up** → **Sign up with GitHub**
3. Authorize Netlify to access your GitHub

#### Step 2: Create New Site

1. Click **Add new site** → **Import an existing project**
2. Choose **GitHub**
3. Find and select `mhd-nanofluid-simulation`

#### Step 3: Configure Build Settings

Netlify should auto-detect these, but verify:

| Setting | Value |
|---------|-------|
| Build command | `npm run build` |
| Publish directory | `build` |

#### Step 4: Deploy

1. Click **Deploy site**
2. Wait 2-5 minutes for the build
3. You'll get a URL like: `https://random-name-123.netlify.app`

#### Step 5: Custom Domain (Optional)

1. In Netlify, go to **Site settings** → **Domain management**
2. Click **Add custom domain** or **Edit site name** to get a better URL

### Updating Your Site

Whenever you make changes:

```bash
git add .
git commit -m "Description of changes"
git push
```

Netlify will automatically rebuild and deploy!

---

## 🎨 Adding Your Own Content

### Adding Videos

1. Create a folder: `public/videos/`
2. Add your video files (MP4 recommended)
3. In `App.js`, find the `renderVideos` function
4. Replace the placeholder with:
   ```jsx
   <video controls width="100%">
     <source src="/videos/your-video-name.mp4" type="video/mp4" />
   </video>
   ```

### Adding Figures/Images

1. Create a folder: `public/images/`
2. Add your image files (PNG or JPG)
3. In `App.js`, find the `renderFigures` function
4. Replace the placeholder with:
   ```jsx
   <img src="/images/Figure_1.png" alt="Figure 1" />
   ```

### Example: Adding Figure 1

In the `renderFigures` function, change:
```jsx
<div className="gallery-item-placeholder">
  <Image />
</div>
```

To:
```jsx
<img 
  src="/images/Figure_1_Grid_Convergence.png" 
  alt="Grid Convergence Study"
  style={{ width: '100%', height: '100%', objectFit: 'cover' }}
/>
```

---

## 📁 Project Structure

```
mhd-nanofluid-app/
├── public/
│   ├── index.html          # Main HTML file
│   ├── images/             # Add your figures here
│   └── videos/             # Add your videos here
├── src/
│   ├── App.js              # Main application component
│   ├── index.js            # Entry point
│   └── index.css           # All styling
├── package.json            # Dependencies and scripts
└── README.md               # This file
```

---

## 📐 Physics Equations

The simulation solves these governing equations:

### Momentum Equation
```
A₁·W'' - A₂·Ha²·W + G = 0
```

### Energy Equation
```
A₃·θ'' + A₁·Pr·Ec·(W')² + A₂·Pr·Ec·Ha²·W² = 0
```

### Boundary Conditions
- **Lower plate (η=0):** W=0, θ=1
- **Upper plate (η=1):** W-λW'=Re, θ'+Bi·θ=0

### Parameters
| Symbol | Name | Description |
|--------|------|-------------|
| Ha | Hartmann Number | Magnetic field strength |
| Re | Reynolds Number | Upper plate velocity |
| Pr | Prandtl Number | Momentum/thermal diffusivity ratio |
| Ec | Eckert Number | Viscous dissipation |
| Bi | Biot Number | Convective heat transfer |
| λ | Slip Parameter | Velocity slip at upper plate |
| G | Pressure Gradient | Axial pressure gradient |
| A₁ | Viscosity Ratio | μnf/μf |
| A₂ | Conductivity Ratio | σnf/σf |
| A₃ | Thermal Ratio | knf/kf |

---

## ❓ Troubleshooting

### "npm is not recognized"
- Node.js is not installed or not in PATH
- Restart your computer after installing Node.js
- Reinstall Node.js from https://nodejs.org

### "Cannot find module..."
- Run `npm install` again
- Delete `node_modules` folder and run `npm install`

### Build Errors on Netlify
- Check the build log for specific errors
- Ensure all files are committed to GitHub
- Verify `package.json` is in the root folder

### App Not Updating After Push
- Wait 2-5 minutes for Netlify to rebuild
- Check Netlify dashboard for build status
- Clear browser cache (Ctrl+Shift+R)

### Graphs Not Showing
- Ensure you have a stable internet connection (charts use recharts library)
- Check browser console for errors (F12 → Console)

---

## 📞 Support

If you encounter issues:
1. Check the troubleshooting section above
2. Search for the error message online
3. Create an issue on the GitHub repository

---

## 📄 License

This project is for educational and research purposes.

---

## 🙏 Acknowledgments

- Nelson Mandela University
- Prof. O.D. Makinde (Supervisor)
- Kigodi et al. (2025) for validation reference

---

**Good luck with your research presentation! 🎓**
