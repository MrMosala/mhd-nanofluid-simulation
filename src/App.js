import React, { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, AreaChart, Area
} from 'recharts';
import { 
  Zap, Thermometer, Play, Image, Video, BookOpen, 
  Droplets, Magnet, Gauge, Activity, Layers, BarChart3, Info
} from 'lucide-react';

// ═══════════════════════════════════════════════════════════════════════════
// PHYSICS ENGINE - MHD NANOFLUID COUETTE FLOW SOLVER
// Based on Proposal_Master.pdf equations
// ═══════════════════════════════════════════════════════════════════════════

/**
 * Solves the coupled MHD Nanofluid Couette Flow equations:
 * Momentum: A1*W'' - A2*Ha²*W + G = 0
 * Energy: A3*θ'' + A1*Pr*Ec*(W')² + A2*Pr*Ec*Ha²*W² = 0
 * 
 * Boundary Conditions:
 * η=0: W=0, θ=1
 * η=1: W-λW'=Re, θ'+Bi*θ=0
 */
function solveMHDCouetteFlow(params) {
  const { A1, A2, A3, Re, Ha, Pr, Ec, Bi, lambda, G, N = 100 } = params;
  
  // Create grid
  const eta = [];
  const h = 1.0 / N;
  for (let i = 0; i <= N; i++) {
    eta.push(i * h);
  }
  
  // Initialize velocity and temperature
  let W = new Array(N + 1).fill(0);
  let Theta = new Array(N + 1).fill(1);
  
  // Set initial guess satisfying BCs
  for (let i = 0; i <= N; i++) {
    W[i] = eta[i] * Re / (1 - lambda);
    Theta[i] = 1 - (Bi / (1 + Bi)) * eta[i];
  }
  
  // Iterative solver (simplified for real-time performance)
  const maxIter = 100;
  const tol = 1e-8;
  
  for (let iter = 0; iter < maxIter; iter++) {
    const W_old = [...W];
    const Theta_old = [...Theta];
    
    // Solve momentum equation using finite differences
    // A1*W'' - A2*Ha²*W + G = 0
    for (let i = 1; i < N; i++) {
      const coeff = A1 / (h * h);
      const diag = 2 * A1 / (h * h) + A2 * Ha * Ha;
      W[i] = (coeff * (W_old[i-1] + W_old[i+1]) + G) / diag;
    }
    
    // Apply velocity BCs
    W[0] = 0; // W(0) = 0
    // W(1) - λW'(1) = Re => W[N] = Re + λ*(W[N] - W[N-1])/h
    W[N] = (Re + lambda * W[N-1] / h) / (1 + lambda / h);
    
    // Calculate velocity derivatives for energy equation
    const Wp = new Array(N + 1).fill(0);
    for (let i = 1; i < N; i++) {
      Wp[i] = (W[i+1] - W[i-1]) / (2 * h);
    }
    Wp[0] = (W[1] - W[0]) / h;
    Wp[N] = (W[N] - W[N-1]) / h;
    
    // Solve energy equation using finite differences
    // A3*θ'' + A1*Pr*Ec*(W')² + A2*Pr*Ec*Ha²*W² = 0
    for (let i = 1; i < N; i++) {
      const source = A1 * Pr * Ec * Wp[i] * Wp[i] + A2 * Pr * Ec * Ha * Ha * W[i] * W[i];
      const coeff = A3 / (h * h);
      const diag = 2 * A3 / (h * h);
      Theta[i] = (coeff * (Theta_old[i-1] + Theta_old[i+1]) + source) / diag;
    }
    
    // Apply temperature BCs
    Theta[0] = 1; // θ(0) = 1
    // θ'(1) + Bi*θ(1) = 0 => (Theta[N] - Theta[N-1])/h + Bi*Theta[N] = 0
    Theta[N] = Theta[N-1] / (1 + h * Bi);
    
    // Check convergence
    let maxDiff = 0;
    for (let i = 0; i <= N; i++) {
      maxDiff = Math.max(maxDiff, Math.abs(W[i] - W_old[i]), Math.abs(Theta[i] - Theta_old[i]));
    }
    if (maxDiff < tol) break;
  }
  
  // Calculate derivatives
  const Wp = new Array(N + 1).fill(0);
  const Thetap = new Array(N + 1).fill(0);
  
  for (let i = 1; i < N; i++) {
    Wp[i] = (W[i+1] - W[i-1]) / (2 * h);
    Thetap[i] = (Theta[i+1] - Theta[i-1]) / (2 * h);
  }
  Wp[0] = (W[1] - W[0]) / h;
  Wp[N] = (W[N] - W[N-1]) / h;
  Thetap[0] = (Theta[1] - Theta[0]) / h;
  Thetap[N] = (Theta[N] - Theta[N-1]) / h;
  
  // Engineering quantities
  const Cf_lower = A1 * Wp[0];
  const Cf_upper = A1 * Wp[N];
  const Nu_lower = -A3 * Thetap[0];
  const Nu_upper = -A3 * Thetap[N];
  
  // Calculate entropy generation
  const Ns = [];
  const Be = [];
  const Ns_heat = [];
  const Ns_fluid = [];
  const Ns_magnetic = [];
  
  for (let i = 0; i <= N; i++) {
    const theta_safe = Math.max(Theta[i], 0.01);
    const ns_h = A3 * (Thetap[i] * Thetap[i]) / (theta_safe * theta_safe);
    const ns_f = A1 * Ec * Pr * (Wp[i] * Wp[i]) / theta_safe;
    const ns_m = A2 * Ec * Pr * Ha * Ha * (W[i] * W[i]) / theta_safe;
    const ns_total = ns_h + ns_f + ns_m;
    
    Ns_heat.push(ns_h);
    Ns_fluid.push(ns_f);
    Ns_magnetic.push(ns_m);
    Ns.push(ns_total);
    Be.push(ns_h / (ns_total + 1e-12));
  }
  
  // Prepare chart data
  const chartData = eta.map((e, i) => ({
    eta: e,
    W: W[i],
    Theta: Theta[i],
    Wp: Wp[i],
    Thetap: Thetap[i],
    Ns: Ns[i],
    Be: Be[i],
    Ns_heat: Ns_heat[i],
    Ns_fluid: Ns_fluid[i],
    Ns_magnetic: Ns_magnetic[i]
  }));
  
  return {
    eta,
    W,
    Theta,
    Wp,
    Thetap,
    Cf_lower,
    Cf_upper,
    Nu_lower,
    Nu_upper,
    Ns,
    Be,
    chartData
  };
}

// ═══════════════════════════════════════════════════════════════════════════
// CUSTOM COMPONENTS
// ═══════════════════════════════════════════════════════════════════════════

// Custom Tooltip for Charts
const CustomTooltip = ({ active, payload, label }) => {
  if (active && payload && payload.length) {
    return (
      <div style={{
        background: 'rgba(10, 14, 23, 0.95)',
        border: '1px solid rgba(0, 212, 255, 0.3)',
        borderRadius: '8px',
        padding: '12px',
        backdropFilter: 'blur(10px)'
      }}>
        <p style={{ color: '#00d4ff', fontFamily: 'Orbitron', fontSize: '0.8rem', marginBottom: '8px' }}>
          η = {typeof label === 'number' ? label.toFixed(3) : label}
        </p>
        {payload.map((item, index) => (
          <p key={index} style={{ color: item.color, fontSize: '0.85rem', margin: '4px 0' }}>
            {item.name}: {typeof item.value === 'number' ? item.value.toFixed(4) : item.value}
          </p>
        ))}
      </div>
    );
  }
  return null;
};

// Parameter Slider Component
const ParameterSlider = ({ label, value, onChange, min, max, step, unit, description }) => {
  return (
    <div className="slider-control">
      <div className="slider-label">
        <span title={description}>{label}</span>
        <span className="slider-value">{value.toFixed(step < 0.01 ? 3 : 2)}{unit}</span>
      </div>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))}
      />
    </div>
  );
};

// Flow Visualization Canvas Component
const FlowVisualization = ({ params, solution }) => {
  const canvasRef = useRef(null);
  const animationRef = useRef(null);
  const particlesRef = useRef([]);
  
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    const width = canvas.width = canvas.offsetWidth * 2;
    const height = canvas.height = canvas.offsetHeight * 2;
    ctx.scale(2, 2);
    
    // Initialize particles
    const numParticles = 150;
    if (particlesRef.current.length === 0) {
      for (let i = 0; i < numParticles; i++) {
        particlesRef.current.push({
          x: Math.random() * (width / 2),
          y: Math.random() * (height / 2),
          size: 1 + Math.random() * 2,
          alpha: 0.3 + Math.random() * 0.5
        });
      }
    }
    
    const animate = () => {
      ctx.fillStyle = 'rgba(10, 14, 23, 0.15)';
      ctx.fillRect(0, 0, width / 2, height / 2);
      
      // Draw plates
      ctx.fillStyle = 'rgba(0, 212, 255, 0.3)';
      ctx.fillRect(0, 0, width / 2, 3);
      ctx.fillStyle = 'rgba(255, 0, 110, 0.3)';
      ctx.fillRect(0, height / 2 - 3, width / 2, 3);
      
      // Draw magnetic field lines
      ctx.strokeStyle = 'rgba(255, 215, 0, 0.1)';
      ctx.lineWidth = 1;
      for (let x = 20; x < width / 2; x += 40) {
        ctx.beginPath();
        ctx.moveTo(x, 10);
        ctx.lineTo(x, height / 2 - 10);
        ctx.stroke();
      }
      
      // Update and draw particles
      particlesRef.current.forEach((p, i) => {
        // Get velocity based on y position
        const eta = 1 - (p.y / (height / 2));
        const etaClamped = Math.max(0, Math.min(1, eta));
        const idx = Math.floor(etaClamped * (solution.W.length - 1));
        const velocity = solution.W[idx] || 0;
        
        // Get temperature for color
        const temp = solution.Theta[idx] || 1;
        
        // Update position
        p.x += velocity * 0.5 + 0.2;
        
        // Wrap around
        if (p.x > width / 2) {
          p.x = 0;
          p.y = Math.random() * (height / 2);
        }
        
        // Draw particle with temperature-based color
        const r = Math.floor(255 * temp);
        const g = Math.floor(100 * (1 - temp) + 150 * temp);
        const b = Math.floor(255 * (1 - temp));
        
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.size, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${p.alpha})`;
        ctx.fill();
        
        // Add glow effect
        const gradient = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, p.size * 3);
        gradient.addColorStop(0, `rgba(${r}, ${g}, ${b}, ${p.alpha * 0.5})`);
        gradient.addColorStop(1, 'transparent');
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.size * 3, 0, Math.PI * 2);
        ctx.fillStyle = gradient;
        ctx.fill();
      });
      
      // Draw labels
      ctx.font = '11px Orbitron';
      ctx.fillStyle = 'rgba(0, 212, 255, 0.8)';
      ctx.fillText('Upper Plate (Moving)', 10, 20);
      ctx.fillStyle = 'rgba(255, 0, 110, 0.8)';
      ctx.fillText('Lower Plate (Stationary)', 10, height / 2 - 10);
      
      // Draw Ha indicator
      ctx.fillStyle = 'rgba(255, 215, 0, 0.8)';
      ctx.fillText(`Ha = ${params.Ha.toFixed(1)}`, width / 2 - 80, 20);
      
      animationRef.current = requestAnimationFrame(animate);
    };
    
    animate();
    
    return () => {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    };
  }, [params, solution]);
  
  return (
    <div className="flow-viz-container">
      <canvas ref={canvasRef} style={{ width: '100%', height: '100%' }} />
    </div>
  );
};

// ═══════════════════════════════════════════════════════════════════════════
// MAIN APP COMPONENT
// ═══════════════════════════════════════════════════════════════════════════

function App() {
  // Navigation State
  const [activeTab, setActiveTab] = useState('simulation');
  
  // Simulation Parameters
  const [params, setParams] = useState({
    // Nanofluid Properties
    A1: 1.2,  // μnf/μf (viscosity ratio)
    A2: 1.5,  // σnf/σf (electrical conductivity ratio)
    A3: 1.3,  // knf/kf (thermal conductivity ratio)
    
    // Flow Parameters
    Re: 1.0,    // Reynolds number
    Ha: 2.0,    // Hartmann number
    G: 0.5,     // Pressure gradient
    lambda: 0.1, // Slip parameter
    
    // Thermal Parameters
    Pr: 6.2,    // Prandtl number
    Ec: 0.01,   // Eckert number
    Bi: 0.5,    // Biot number
    
    // Numerical
    N: 100      // Grid points
  });
  
  // Solution State
  const solution = useMemo(() => {
    return solveMHDCouetteFlow(params);
  }, [params]);
  
  // Parameter Update Handler
  const updateParam = useCallback((key, value) => {
    setParams(prev => ({ ...prev, [key]: value }));
  }, []);
  
  // Tab Content Renderers
  const renderSimulation = () => (
    <div className="simulation-grid">
      {/* Controls Panel */}
      <div className="controls-panel">
        <div className="card">
          <div className="card-header">
            <div className="icon"><Droplets size={20} /></div>
            <h2>Nanofluid Properties</h2>
          </div>
          <div className="param-group">
            <div className="param-group-title">Volume Fraction Effects</div>
            <ParameterSlider
              label="A₁ (μnf/μf)"
              value={params.A1}
              onChange={(v) => updateParam('A1', v)}
              min={1.0} max={2.0} step={0.05}
              unit=""
              description="Viscosity ratio"
            />
            <ParameterSlider
              label="A₂ (σnf/σf)"
              value={params.A2}
              onChange={(v) => updateParam('A2', v)}
              min={1.0} max={3.0} step={0.05}
              unit=""
              description="Electrical conductivity ratio"
            />
            <ParameterSlider
              label="A₃ (knf/kf)"
              value={params.A3}
              onChange={(v) => updateParam('A3', v)}
              min={1.0} max={2.0} step={0.05}
              unit=""
              description="Thermal conductivity ratio"
            />
          </div>
        </div>
        
        <div className="card">
          <div className="card-header">
            <div className="icon"><Magnet size={20} /></div>
            <h2>MHD Parameters</h2>
          </div>
          <div className="param-group">
            <div className="param-group-title">Magnetic & Flow</div>
            <ParameterSlider
              label="Ha (Hartmann)"
              value={params.Ha}
              onChange={(v) => updateParam('Ha', v)}
              min={0} max={10} step={0.1}
              unit=""
              description="Magnetic field strength"
            />
            <ParameterSlider
              label="Re (Reynolds)"
              value={params.Re}
              onChange={(v) => updateParam('Re', v)}
              min={0.1} max={5} step={0.1}
              unit=""
              description="Upper plate velocity"
            />
            <ParameterSlider
              label="G (Pressure)"
              value={params.G}
              onChange={(v) => updateParam('G', v)}
              min={0} max={2} step={0.1}
              unit=""
              description="Pressure gradient parameter"
            />
            <ParameterSlider
              label="λ (Slip)"
              value={params.lambda}
              onChange={(v) => updateParam('lambda', v)}
              min={0} max={0.5} step={0.01}
              unit=""
              description="Velocity slip parameter"
            />
          </div>
        </div>
        
        <div className="card">
          <div className="card-header">
            <div className="icon"><Thermometer size={20} /></div>
            <h2>Thermal Parameters</h2>
          </div>
          <div className="param-group">
            <div className="param-group-title">Heat Transfer</div>
            <ParameterSlider
              label="Pr (Prandtl)"
              value={params.Pr}
              onChange={(v) => updateParam('Pr', v)}
              min={0.7} max={20} step={0.1}
              unit=""
              description="Momentum to thermal diffusivity ratio"
            />
            <ParameterSlider
              label="Ec (Eckert)"
              value={params.Ec}
              onChange={(v) => updateParam('Ec', v)}
              min={0} max={0.2} step={0.005}
              unit=""
              description="Viscous dissipation"
            />
            <ParameterSlider
              label="Bi (Biot)"
              value={params.Bi}
              onChange={(v) => updateParam('Bi', v)}
              min={0.1} max={5} step={0.1}
              unit=""
              description="Convective heat transfer coefficient"
            />
          </div>
        </div>
      </div>
      
      {/* Visualization Panel */}
      <div className="visualization-panel">
        {/* Flow Animation */}
        <FlowVisualization params={params} solution={solution} />
        
        {/* Results Summary */}
        <div className="results-panel">
          <div className="result-card cyan">
            <div className="label">Cf (Lower)</div>
            <div className="value">{solution.Cf_lower.toFixed(4)}</div>
          </div>
          <div className="result-card magenta">
            <div className="label">Cf (Upper)</div>
            <div className="value">{solution.Cf_upper.toFixed(4)}</div>
          </div>
          <div className="result-card gold">
            <div className="label">Nu (Lower)</div>
            <div className="value">{solution.Nu_lower.toFixed(4)}</div>
          </div>
          <div className="result-card emerald">
            <div className="label">Nu (Upper)</div>
            <div className="value">{solution.Nu_upper.toFixed(4)}</div>
          </div>
        </div>
        
        {/* Charts */}
        <div className="charts-grid">
          {/* Velocity Profile */}
          <div className="chart-container">
            <div className="chart-title">
              <span className="dot"></span>
              Velocity Profile W(η)
            </div>
            <ResponsiveContainer width="100%" height={220}>
              <LineChart data={solution.chartData} margin={{ top: 10, right: 30, left: 0, bottom: 0 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis 
                  dataKey="eta" 
                  stroke="rgba(255,255,255,0.5)"
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }}
                  label={{ value: 'η', position: 'insideBottomRight', offset: -5, fill: '#00d4ff' }}
                />
                <YAxis 
                  stroke="rgba(255,255,255,0.5)"
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }}
                />
                <Tooltip content={<CustomTooltip />} />
                <Line 
                  type="monotone" 
                  dataKey="W" 
                  stroke="#00d4ff" 
                  strokeWidth={2}
                  dot={false}
                  name="Velocity W"
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
          
          {/* Temperature Profile */}
          <div className="chart-container">
            <div className="chart-title">
              <span className="dot magenta"></span>
              Temperature Profile θ(η)
            </div>
            <ResponsiveContainer width="100%" height={220}>
              <LineChart data={solution.chartData} margin={{ top: 10, right: 30, left: 0, bottom: 0 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis 
                  dataKey="eta" 
                  stroke="rgba(255,255,255,0.5)"
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }}
                  label={{ value: 'η', position: 'insideBottomRight', offset: -5, fill: '#ff006e' }}
                />
                <YAxis 
                  stroke="rgba(255,255,255,0.5)"
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }}
                />
                <Tooltip content={<CustomTooltip />} />
                <Line 
                  type="monotone" 
                  dataKey="Theta" 
                  stroke="#ff006e" 
                  strokeWidth={2}
                  dot={false}
                  name="Temperature θ"
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
          
          {/* Entropy Generation */}
          <div className="chart-container">
            <div className="chart-title">
              <span className="dot gold"></span>
              Entropy Generation Ns(η)
            </div>
            <ResponsiveContainer width="100%" height={220}>
              <AreaChart data={solution.chartData} margin={{ top: 10, right: 30, left: 0, bottom: 0 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis 
                  dataKey="eta" 
                  stroke="rgba(255,255,255,0.5)"
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }}
                />
                <YAxis 
                  stroke="rgba(255,255,255,0.5)"
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }}
                />
                <Tooltip content={<CustomTooltip />} />
                <Area 
                  type="monotone" 
                  dataKey="Ns_heat" 
                  stackId="1"
                  stroke="#ff006e" 
                  fill="rgba(255, 0, 110, 0.3)"
                  name="Heat Transfer"
                />
                <Area 
                  type="monotone" 
                  dataKey="Ns_fluid" 
                  stackId="1"
                  stroke="#00d4ff" 
                  fill="rgba(0, 212, 255, 0.3)"
                  name="Fluid Friction"
                />
                <Area 
                  type="monotone" 
                  dataKey="Ns_magnetic" 
                  stackId="1"
                  stroke="#ffd700" 
                  fill="rgba(255, 215, 0, 0.3)"
                  name="Magnetic"
                />
              </AreaChart>
            </ResponsiveContainer>
          </div>
          
          {/* Bejan Number */}
          <div className="chart-container">
            <div className="chart-title">
              <span className="dot emerald"></span>
              Bejan Number Be(η)
            </div>
            <ResponsiveContainer width="100%" height={220}>
              <LineChart data={solution.chartData} margin={{ top: 10, right: 30, left: 0, bottom: 0 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis 
                  dataKey="eta" 
                  stroke="rgba(255,255,255,0.5)"
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }}
                />
                <YAxis 
                  stroke="rgba(255,255,255,0.5)"
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }}
                  domain={[0, 1]}
                />
                <Tooltip content={<CustomTooltip />} />
                <Line 
                  type="monotone" 
                  dataKey="Be" 
                  stroke="#00ff9f" 
                  strokeWidth={2}
                  dot={false}
                  name="Bejan Number"
                />
                {/* Reference line at Be = 0.5 */}
                <Line 
                  type="monotone" 
                  data={[{eta: 0, ref: 0.5}, {eta: 1, ref: 0.5}]}
                  dataKey="ref"
                  stroke="rgba(255,255,255,0.3)"
                  strokeDasharray="5 5"
                  strokeWidth={1}
                  dot={false}
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </div>
      </div>
    </div>
  );
  
  const renderVideos = () => (
    <div className="gallery-grid">
      {[
        { title: "Flow Animation with Ha Variation", desc: "Watch how magnetic field strength affects velocity profiles" },
        { title: "Temperature Field Evolution", desc: "Visualization of thermal boundary layer development" },
        { title: "Nanoparticle Concentration Effects", desc: "Impact of volume fraction on heat transfer" },
        { title: "Entropy Generation Analysis", desc: "Thermodynamic irreversibility breakdown" },
        { title: "Parametric Study Overview", desc: "Comprehensive parameter sensitivity analysis" },
        { title: "SQLM Numerical Method", desc: "Spectral quasi-linearization method explained" }
      ].map((video, index) => (
        <div key={index} className="video-card">
          <div className="video-placeholder">
            <Video />
            <p>Video placeholder - Upload your video</p>
          </div>
          <div className="video-info">
            <h3>{video.title}</h3>
            <p>{video.desc}</p>
            <p style={{ marginTop: '1rem', color: 'var(--accent-cyan)', fontSize: '0.85rem' }}>
              <Info size={14} style={{ display: 'inline', marginRight: '4px' }} />
              Replace this placeholder with your actual video file
            </p>
          </div>
        </div>
      ))}
    </div>
  );
  
  const renderFigures = () => (
    <div className="gallery-grid">
      {[
        { title: "Grid Convergence Study", desc: "Spectral accuracy validation showing exponential convergence" },
        { title: "Analytical Validation", desc: "Comparison with exact solutions for limiting cases" },
        { title: "Hartmann Number Effects", desc: "Velocity and temperature profiles for varying Ha" },
        { title: "Reynolds Number Study", desc: "Flow characteristics at different Re values" },
        { title: "Prandtl Number Analysis", desc: "Thermal boundary layer response to Pr variations" },
        { title: "Eckert Number Effects", desc: "Viscous dissipation impact on temperature" },
        { title: "3D Ha-Re Surface", desc: "Skin friction coefficient parameter space" },
        { title: "3D Pr-Ec Surface", desc: "Nusselt number thermal parameter coupling" },
        { title: "Entropy Distribution", desc: "Spatial distribution of irreversibility sources" },
        { title: "Bejan Number Profiles", desc: "Heat transfer vs friction dominance analysis" },
        { title: "Nanofluid Enhancement", desc: "Comparison of nanofluid vs base fluid performance" },
        { title: "Literature Validation", desc: "Comparison with Kigodi et al. (2025) results" }
      ].map((figure, index) => (
        <div key={index} className="gallery-item">
          <div className="gallery-item-image">
            <div className="gallery-item-placeholder">
              <Image />
            </div>
          </div>
          <div className="gallery-item-info">
            <h3>Figure {index + 1}: {figure.title}</h3>
            <p>{figure.desc}</p>
            <p style={{ marginTop: '0.5rem', color: 'var(--accent-cyan)', fontSize: '0.8rem' }}>
              Replace with: Figure_{index + 1}.png
            </p>
          </div>
        </div>
      ))}
    </div>
  );
  
  const renderTheory = () => (
    <div className="theory-content">
      <div className="equation-card">
        <h3><Activity size={20} /> Governing Equations</h3>
        <div className="equation">
          <strong>Momentum Equation:</strong><br/>
          A₁·W'' - A₂·Ha²·W + G = 0
        </div>
        <p className="equation-description">
          Where W is the dimensionless velocity, Ha is the Hartmann number representing 
          magnetic field strength, and G is the pressure gradient parameter. The A₁ and A₂ 
          coefficients account for nanofluid viscosity and electrical conductivity enhancement.
        </p>
        
        <div className="equation">
          <strong>Energy Equation:</strong><br/>
          A₃·θ'' + A₁·Pr·Ec·(W')² + A₂·Pr·Ec·Ha²·W² = 0
        </div>
        <p className="equation-description">
          The energy equation includes viscous dissipation (Eckert number Ec) and Joule 
          heating (Ha² term). The Prandtl number Pr relates momentum and thermal diffusivity.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><Layers size={20} /> Boundary Conditions</h3>
        <div className="equation">
          <strong>At η = 0 (Lower Plate):</strong><br/>
          W(0) = 0, θ(0) = 1
        </div>
        <p className="equation-description">
          The lower plate is stationary with a fixed wall temperature.
        </p>
        
        <div className="equation">
          <strong>At η = 1 (Upper Plate):</strong><br/>
          W(1) - λ·W'(1) = Re<br/>
          θ'(1) + Bi·θ(1) = 0
        </div>
        <p className="equation-description">
          The upper plate moves with velocity Re, with slip parameter λ. 
          Convective cooling occurs with Biot number Bi.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><Gauge size={20} /> Engineering Quantities</h3>
        <div className="equation">
          <strong>Skin Friction:</strong><br/>
          Cf = A₁ · dW/dη |<sub>η=0,1</sub>
        </div>
        <div className="equation">
          <strong>Nusselt Number:</strong><br/>
          Nu = -A₃ · dθ/dη |<sub>η=0,1</sub>
        </div>
        <p className="equation-description">
          These quantities measure wall shear stress and heat transfer rate at both plates.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><BarChart3 size={20} /> Entropy Generation</h3>
        <div className="equation">
          <strong>Total Entropy:</strong><br/>
          Ns = Ns,heat + Ns,fluid + Ns,magnetic
        </div>
        <div className="equation">
          <strong>Bejan Number:</strong><br/>
          Be = Ns,heat / Ns,total
        </div>
        <p className="equation-description">
          Entropy generation analysis reveals thermodynamic irreversibility sources. 
          Be {'>'} 0.5 indicates heat transfer dominance; Be {'<'} 0.5 indicates friction/magnetic dominance.
        </p>
      </div>
      
      <div className="equation-card" style={{ gridColumn: 'span 2' }}>
        <h3><Droplets size={20} /> Nanofluid Properties</h3>
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem' }}>
          <div className="equation">
            <strong>Density:</strong><br/>
            ρnf = (1-φ)ρf + φρs
          </div>
          <div className="equation">
            <strong>Heat Capacity:</strong><br/>
            (ρCp)nf = (1-φ)(ρCp)f + φ(ρCp)s
          </div>
          <div className="equation">
            <strong>Viscosity (Brinkman):</strong><br/>
            μnf = μf / (1-φ)^2.5
          </div>
          <div className="equation">
            <strong>Thermal Conductivity (Maxwell):</strong><br/>
            knf/kf = (ks+2kf-2φ(kf-ks)) / (ks+2kf+φ(kf-ks))
          </div>
        </div>
        <p className="equation-description" style={{ marginTop: '1rem' }}>
          These correlations model the effective thermophysical properties of the nanofluid 
          as functions of nanoparticle volume fraction φ. The A₁, A₂, A₃ parameters in the 
          governing equations represent ratios of nanofluid to base fluid properties.
        </p>
      </div>
    </div>
  );
  
  return (
    <div className="app">
      {/* Header */}
      <header className="header">
        <div className="header-content">
          <div className="logo-section">
            <div className="logo-icon">
              <Zap size={28} />
            </div>
            <div className="logo-text">
              <h1>MHD Nanofluid Flow</h1>
              <p>Couette Flow Simulation</p>
            </div>
          </div>
          
          <nav className="nav-tabs">
            <button 
              className={`nav-tab ${activeTab === 'simulation' ? 'active' : ''}`}
              onClick={() => setActiveTab('simulation')}
            >
              <Play size={16} /> Simulation
            </button>
            <button 
              className={`nav-tab ${activeTab === 'videos' ? 'active' : ''}`}
              onClick={() => setActiveTab('videos')}
            >
              <Video size={16} /> Videos
            </button>
            <button 
              className={`nav-tab ${activeTab === 'figures' ? 'active' : ''}`}
              onClick={() => setActiveTab('figures')}
            >
              <Image size={16} /> Figures
            </button>
            <button 
              className={`nav-tab ${activeTab === 'theory' ? 'active' : ''}`}
              onClick={() => setActiveTab('theory')}
            >
              <BookOpen size={16} /> Theory
            </button>
          </nav>
        </div>
      </header>
      
      {/* Main Content */}
      <main className="main-content">
        {activeTab === 'simulation' && renderSimulation()}
        {activeTab === 'videos' && renderVideos()}
        {activeTab === 'figures' && renderFigures()}
        {activeTab === 'theory' && renderTheory()}
      </main>
      
      {/* Footer */}
      <footer className="footer">
        <p>
          <strong>Research Project:</strong> Thermal and Magnetohydrodynamic Analysis of Nanofluid Couette Flow
          <br />
          <strong>Candidate:</strong> Mr. S.I. Mosala | <strong>Supervisor:</strong> Prof. O.D. Makinde
          <br />
          Nelson Mandela University | December 2025
        </p>
      </footer>
    </div>
  );
}

export default App;
