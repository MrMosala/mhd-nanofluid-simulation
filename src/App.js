import React, { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, AreaChart, Area
} from 'recharts';
import { 
  Zap, Thermometer, Play, Image, Video, BookOpen, 
  Droplets, Magnet, Gauge, Activity, Layers, BarChart3, Info,
  Sliders, X, ChevronDown, Wind, TrendingUp
} from 'lucide-react';

// ═══════════════════════════════════════════════════════════════════════════
// PHYSICS ENGINE - MHD NANOFLUID COUETTE FLOW SOLVER
// Based on Proposal_Master.pdf equations
// ═══════════════════════════════════════════════════════════════════════════

function solveMHDCouetteFlow(params) {
  const { A1, A2, A3, Re, Ha, Pr, Ec, Bi, lambda, G, N = 100 } = params;
  
  const eta = [];
  const h = 1.0 / N;
  for (let i = 0; i <= N; i++) {
    eta.push(i * h);
  }
  
  let W = new Array(N + 1).fill(0);
  let Theta = new Array(N + 1).fill(1);
  
  for (let i = 0; i <= N; i++) {
    W[i] = eta[i] * Re / (1 - lambda);
    Theta[i] = 1 - (Bi / (1 + Bi)) * eta[i];
  }
  
  const maxIter = 100;
  const tol = 1e-8;
  
  for (let iter = 0; iter < maxIter; iter++) {
    const W_old = [...W];
    const Theta_old = [...Theta];
    
    for (let i = 1; i < N; i++) {
      const coeff = A1 / (h * h);
      const diag = 2 * A1 / (h * h) + A2 * Ha * Ha;
      W[i] = (coeff * (W_old[i-1] + W_old[i+1]) + G) / diag;
    }
    
    W[0] = 0;
    W[N] = (Re + lambda * W[N-1] / h) / (1 + lambda / h);
    
    const Wp = new Array(N + 1).fill(0);
    for (let i = 1; i < N; i++) {
      Wp[i] = (W[i+1] - W[i-1]) / (2 * h);
    }
    Wp[0] = (W[1] - W[0]) / h;
    Wp[N] = (W[N] - W[N-1]) / h;
    
    for (let i = 1; i < N; i++) {
      const source = A1 * Pr * Ec * Wp[i] * Wp[i] + A2 * Pr * Ec * Ha * Ha * W[i] * W[i];
      const coeff = A3 / (h * h);
      const diag = 2 * A3 / (h * h);
      Theta[i] = (coeff * (Theta_old[i-1] + Theta_old[i+1]) + source) / diag;
    }
    
    Theta[0] = 1;
    Theta[N] = Theta[N-1] / (1 + h * Bi);
    
    let maxDiff = 0;
    for (let i = 0; i <= N; i++) {
      maxDiff = Math.max(maxDiff, Math.abs(W[i] - W_old[i]), Math.abs(Theta[i] - Theta_old[i]));
    }
    if (maxDiff < tol) break;
  }
  
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
  
  const Cf_lower = A1 * Wp[0];
  const Cf_upper = A1 * Wp[N];
  const Nu_lower = -A3 * Thetap[0];
  const Nu_upper = -A3 * Thetap[N];
  
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
    eta, W, Theta, Wp, Thetap,
    Cf_lower, Cf_upper, Nu_lower, Nu_upper,
    Ns, Be, chartData,
    maxW: Math.max(...W),
    minTheta: Math.min(...Theta),
    maxTheta: Math.max(...Theta),
    avgBe: Be.reduce((a, b) => a + b, 0) / Be.length
  };
}

// ═══════════════════════════════════════════════════════════════════════════
// CUSTOM COMPONENTS
// ═══════════════════════════════════════════════════════════════════════════

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

// Parameter Accordion Component
const ParamAccordion = ({ title, icon: Icon, children, defaultOpen = false }) => {
  const [isOpen, setIsOpen] = useState(defaultOpen);
  
  return (
    <div className="param-accordion">
      <div 
        className={`param-accordion-header ${isOpen ? 'open' : ''}`}
        onClick={() => setIsOpen(!isOpen)}
      >
        <div className="param-accordion-title">
          <Icon size={18} />
          {title}
        </div>
        <ChevronDown className={`param-accordion-arrow ${isOpen ? 'open' : ''}`} size={18} />
      </div>
      <div className={`param-accordion-content ${isOpen ? 'open' : ''}`}>
        <div className="param-accordion-inner">
          {children}
        </div>
      </div>
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
    const rect = canvas.getBoundingClientRect();
    const width = canvas.width = rect.width * 2;
    const height = canvas.height = rect.height * 2;
    ctx.scale(2, 2);
    
    const numParticles = 120;
    particlesRef.current = [];
    for (let i = 0; i < numParticles; i++) {
      particlesRef.current.push({
        x: Math.random() * (width / 2),
        y: Math.random() * (height / 2),
        size: 1.5 + Math.random() * 2.5,
        alpha: 0.4 + Math.random() * 0.5
      });
    }
    
    const animate = () => {
      ctx.fillStyle = 'rgba(10, 14, 23, 0.12)';
      ctx.fillRect(0, 0, width / 2, height / 2);
      
      // Draw plates
      const gradient1 = ctx.createLinearGradient(0, 0, width / 2, 0);
      gradient1.addColorStop(0, 'rgba(0, 212, 255, 0.6)');
      gradient1.addColorStop(1, 'rgba(0, 212, 255, 0.2)');
      ctx.fillStyle = gradient1;
      ctx.fillRect(0, 0, width / 2, 4);
      
      const gradient2 = ctx.createLinearGradient(0, 0, width / 2, 0);
      gradient2.addColorStop(0, 'rgba(255, 0, 110, 0.6)');
      gradient2.addColorStop(1, 'rgba(255, 0, 110, 0.2)');
      ctx.fillStyle = gradient2;
      ctx.fillRect(0, height / 2 - 4, width / 2, 4);
      
      // Draw magnetic field lines
      ctx.strokeStyle = 'rgba(255, 215, 0, 0.08)';
      ctx.lineWidth = 1;
      ctx.setLineDash([5, 10]);
      for (let x = 30; x < width / 2; x += 50) {
        ctx.beginPath();
        ctx.moveTo(x, 15);
        ctx.lineTo(x, height / 2 - 15);
        ctx.stroke();
      }
      ctx.setLineDash([]);
      
      // Update and draw particles
      particlesRef.current.forEach((p) => {
        const eta = 1 - (p.y / (height / 2));
        const etaClamped = Math.max(0, Math.min(1, eta));
        const idx = Math.floor(etaClamped * (solution.W.length - 1));
        const velocity = solution.W[idx] || 0;
        const temp = solution.Theta[idx] || 1;
        
        p.x += velocity * 0.6 + 0.3;
        
        if (p.x > width / 2) {
          p.x = 0;
          p.y = Math.random() * (height / 2);
        }
        
        const r = Math.floor(255 * temp);
        const g = Math.floor(100 + 150 * (1 - temp));
        const b = Math.floor(255 * (1 - temp * 0.7));
        
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.size, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${p.alpha})`;
        ctx.fill();
        
        const glowGradient = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, p.size * 4);
        glowGradient.addColorStop(0, `rgba(${r}, ${g}, ${b}, ${p.alpha * 0.4})`);
        glowGradient.addColorStop(1, 'transparent');
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.size * 4, 0, Math.PI * 2);
        ctx.fillStyle = glowGradient;
        ctx.fill();
      });
      
      // Labels
      ctx.font = '12px Orbitron';
      ctx.fillStyle = 'rgba(0, 212, 255, 0.9)';
      ctx.fillText('Upper Plate (Moving) → Re = ' + params.Re.toFixed(1), 10, 22);
      
      ctx.fillStyle = 'rgba(255, 0, 110, 0.9)';
      ctx.fillText('Lower Plate (Stationary)', 10, height / 2 - 12);
      
      ctx.fillStyle = 'rgba(255, 215, 0, 0.9)';
      ctx.textAlign = 'right';
      ctx.fillText('Ha = ' + params.Ha.toFixed(1), width / 2 - 10, 22);
      ctx.textAlign = 'left';
      
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
  const [activeTab, setActiveTab] = useState('simulation');
  const [controlsOpen, setControlsOpen] = useState(false);
  
  const [params, setParams] = useState({
    A1: 1.2, A2: 1.5, A3: 1.3,
    Re: 1.0, Ha: 2.0, G: 0.5, lambda: 0.1,
    Pr: 6.2, Ec: 0.01, Bi: 0.5,
    N: 100
  });
  
  const solution = useMemo(() => solveMHDCouetteFlow(params), [params]);
  
  const updateParam = useCallback((key, value) => {
    setParams(prev => ({ ...prev, [key]: value }));
  }, []);

  // Floating Controls Panel
  const FloatingControls = () => (
    <>
      <div className="floating-controls">
        <button 
          className="controls-toggle-btn"
          onClick={() => setControlsOpen(true)}
          aria-label="Open Parameters"
        >
          <Sliders />
        </button>
      </div>
      
      <div 
        className={`controls-panel-overlay ${controlsOpen ? 'open' : ''}`}
        onClick={() => setControlsOpen(false)}
      />
      
      <div className={`controls-panel-drawer ${controlsOpen ? 'open' : ''}`}>
        <div className="drawer-handle" />
        <div className="drawer-header">
          <h3>⚙️ Simulation Parameters</h3>
          <button className="drawer-close-btn" onClick={() => setControlsOpen(false)}>
            <X size={18} />
          </button>
        </div>
        <div className="drawer-content">
          <ParamAccordion title="Nanofluid Properties" icon={Droplets} defaultOpen={true}>
            <ParameterSlider label="A₁ (μnf/μf)" value={params.A1} onChange={(v) => updateParam('A1', v)} min={1.0} max={2.0} step={0.05} unit="" />
            <ParameterSlider label="A₂ (σnf/σf)" value={params.A2} onChange={(v) => updateParam('A2', v)} min={1.0} max={3.0} step={0.05} unit="" />
            <ParameterSlider label="A₃ (knf/kf)" value={params.A3} onChange={(v) => updateParam('A3', v)} min={1.0} max={2.0} step={0.05} unit="" />
          </ParamAccordion>
          
          <ParamAccordion title="MHD Parameters" icon={Magnet} defaultOpen={true}>
            <ParameterSlider label="Ha (Hartmann)" value={params.Ha} onChange={(v) => updateParam('Ha', v)} min={0} max={10} step={0.1} unit="" />
            <ParameterSlider label="Re (Reynolds)" value={params.Re} onChange={(v) => updateParam('Re', v)} min={0.1} max={5} step={0.1} unit="" />
            <ParameterSlider label="G (Pressure)" value={params.G} onChange={(v) => updateParam('G', v)} min={0} max={2} step={0.1} unit="" />
            <ParameterSlider label="λ (Slip)" value={params.lambda} onChange={(v) => updateParam('lambda', v)} min={0} max={0.5} step={0.01} unit="" />
          </ParamAccordion>
          
          <ParamAccordion title="Thermal Parameters" icon={Thermometer}>
            <ParameterSlider label="Pr (Prandtl)" value={params.Pr} onChange={(v) => updateParam('Pr', v)} min={0.7} max={20} step={0.1} unit="" />
            <ParameterSlider label="Ec (Eckert)" value={params.Ec} onChange={(v) => updateParam('Ec', v)} min={0} max={0.2} step={0.005} unit="" />
            <ParameterSlider label="Bi (Biot)" value={params.Bi} onChange={(v) => updateParam('Bi', v)} min={0.1} max={5} step={0.1} unit="" />
          </ParamAccordion>
        </div>
      </div>
    </>
  );

  // Results Panel Component
  const ResultsPanel = () => (
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
  );

  // Tab Renderers
  const renderSimulation = () => (
    <div className="visualization-section animate-slide-up">
      <ResultsPanel />
      <FlowVisualization params={params} solution={solution} />
      
      <div className="physics-box">
        <h4><Activity size={18} /> Understanding the Flow Visualization</h4>
        <p>
          This animation shows nanofluid particles flowing between two parallel plates. The 
          <strong style={{color: 'var(--accent-cyan)'}}> upper plate (cyan)</strong> moves with velocity Re, 
          while the <strong style={{color: 'var(--accent-magenta)'}}> lower plate (magenta)</strong> remains stationary.
        </p>
        
        <div className="physics-highlight">
          <strong>Particle Color:</strong> Indicates temperature - 
          <span style={{color: '#ff6b6b'}}> warmer (red)</span> to 
          <span style={{color: '#4dabf7'}}> cooler (blue)</span>
        </div>
        
        <div className="physics-highlight gold">
          <strong>Yellow dashed lines:</strong> Represent the transverse magnetic field (B₀). 
          Increasing Ha (Hartmann number) strengthens the Lorentz force, which opposes fluid motion.
        </div>
        
        <p><strong>Try this:</strong> Increase Ha to see particles slow down due to magnetic damping!</p>
      </div>
    </div>
  );

  const renderVelocity = () => (
    <div className="visualization-section animate-slide-up">
      <ResultsPanel />
      
      <div className="chart-section">
        <div className="chart-header">
          <span className="dot"></span>
          <h3>Velocity Profile W(η)</h3>
        </div>
        <div className="chart-wrapper">
          <ResponsiveContainer width="100%" height="100%">
            <LineChart data={solution.chartData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
              <XAxis 
                dataKey="eta" 
                stroke="rgba(255,255,255,0.5)"
                tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 12 }}
                label={{ value: 'η (dimensionless position)', position: 'insideBottom', offset: -10, fill: '#00d4ff' }}
              />
              <YAxis 
                stroke="rgba(255,255,255,0.5)"
                tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 12 }}
                label={{ value: 'W (velocity)', angle: -90, position: 'insideLeft', fill: '#00d4ff' }}
              />
              <Tooltip content={<CustomTooltip />} />
              <Line type="monotone" dataKey="W" stroke="#00d4ff" strokeWidth={3} dot={false} name="Velocity W" />
            </LineChart>
          </ResponsiveContainer>
        </div>
      </div>
      
      <div className="dual-chart-grid">
        <div className="chart-section">
          <div className="chart-header">
            <span className="dot emerald"></span>
            <h3>Velocity Gradient W'(η)</h3>
          </div>
          <div className="chart-wrapper">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={solution.chartData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis dataKey="eta" stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }} />
                <YAxis stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }} />
                <Tooltip content={<CustomTooltip />} />
                <Line type="monotone" dataKey="Wp" stroke="#00ff9f" strokeWidth={2} dot={false} name="Velocity Gradient W'" />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </div>
        
        <div className="physics-box" style={{ margin: 0 }}>
          <h4><TrendingUp size={18} /> Key Metrics</h4>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem', marginBottom: '1rem' }}>
            <div className="result-card cyan">
              <div className="label">Max Velocity</div>
              <div className="value">{solution.maxW.toFixed(4)}</div>
            </div>
            <div className="result-card emerald">
              <div className="label">W'(0)</div>
              <div className="value">{solution.Wp[0].toFixed(4)}</div>
            </div>
          </div>
          
          <h4><Info size={18} /> Skin Friction Coefficient</h4>
          <p>The skin friction is calculated as:</p>
          <div className="equation-inline">Cf = A₁ × dW/dη</div>
          <p style={{ marginTop: '0.5rem' }}>
            At the lower plate: <strong style={{color: 'var(--accent-cyan)'}}>Cf = {solution.Cf_lower.toFixed(4)}</strong>
          </p>
        </div>
      </div>
      
      <div className="physics-box">
        <h4><Wind size={18} /> Physics of Velocity Distribution</h4>
        
        <div className="physics-highlight">
          <strong>Governing Equation:</strong>
          <div className="equation-inline" style={{ display: 'block', marginTop: '0.5rem' }}>
            A₁·W'' - A₂·Ha²·W + G = 0
          </div>
        </div>
        
        <p><strong>Physical Interpretation:</strong></p>
        <ul>
          <li><strong>A₁·W''</strong> — Viscous diffusion (nanofluid viscosity effect)</li>
          <li><strong>A₂·Ha²·W</strong> — Lorentz force (magnetic damping)</li>
          <li><strong>G</strong> — Pressure gradient driving force</li>
        </ul>
        
        <div className="physics-highlight magenta">
          <strong>Effect of Hartmann Number (Ha):</strong><br/>
          Increasing Ha → Stronger magnetic field → Greater Lorentz force → 
          <em> Velocity decreases</em> as the magnetic field opposes fluid motion.
        </div>
        
        <p><strong>Boundary Conditions:</strong></p>
        <ul>
          <li><strong>η = 0 (Lower plate):</strong> W = 0 (no-slip, stationary)</li>
          <li><strong>η = 1 (Upper plate):</strong> W - λW' = Re (slip condition)</li>
        </ul>
      </div>
    </div>
  );

  const renderTemperature = () => (
    <div className="visualization-section animate-slide-up">
      <ResultsPanel />
      
      <div className="chart-section">
        <div className="chart-header">
          <span className="dot magenta"></span>
          <h3>Temperature Profile θ(η)</h3>
        </div>
        <div className="chart-wrapper">
          <ResponsiveContainer width="100%" height="100%">
            <AreaChart data={solution.chartData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
              <defs>
                <linearGradient id="tempGradient" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#ff006e" stopOpacity={0.4}/>
                  <stop offset="95%" stopColor="#ff006e" stopOpacity={0.05}/>
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
              <XAxis 
                dataKey="eta" 
                stroke="rgba(255,255,255,0.5)"
                tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 12 }}
                label={{ value: 'η (dimensionless position)', position: 'insideBottom', offset: -10, fill: '#ff006e' }}
              />
              <YAxis 
                stroke="rgba(255,255,255,0.5)"
                tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 12 }}
                label={{ value: 'θ (temperature)', angle: -90, position: 'insideLeft', fill: '#ff006e' }}
              />
              <Tooltip content={<CustomTooltip />} />
              <Area type="monotone" dataKey="Theta" stroke="#ff006e" strokeWidth={3} fill="url(#tempGradient)" name="Temperature θ" />
            </AreaChart>
          </ResponsiveContainer>
        </div>
      </div>
      
      <div className="dual-chart-grid">
        <div className="chart-section">
          <div className="chart-header">
            <span className="dot gold"></span>
            <h3>Temperature Gradient θ'(η)</h3>
          </div>
          <div className="chart-wrapper">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={solution.chartData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis dataKey="eta" stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }} />
                <YAxis stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 11 }} />
                <Tooltip content={<CustomTooltip />} />
                <Line type="monotone" dataKey="Thetap" stroke="#ffd700" strokeWidth={2} dot={false} name="Temperature Gradient θ'" />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </div>
        
        <div className="physics-box" style={{ margin: 0 }}>
          <h4><Thermometer size={18} /> Temperature Statistics</h4>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem', marginBottom: '1rem' }}>
            <div className="result-card magenta">
              <div className="label">Max θ</div>
              <div className="value">{solution.maxTheta.toFixed(4)}</div>
            </div>
            <div className="result-card gold">
              <div className="label">Min θ</div>
              <div className="value">{solution.minTheta.toFixed(4)}</div>
            </div>
          </div>
          
          <h4><Info size={18} /> Nusselt Number</h4>
          <p>Heat transfer rate:</p>
          <div className="equation-inline">Nu = -A₃ × dθ/dη</div>
          <p style={{ marginTop: '0.5rem' }}>
            At the lower plate: <strong style={{color: 'var(--accent-gold)'}}>Nu = {solution.Nu_lower.toFixed(4)}</strong>
          </p>
        </div>
      </div>
      
      <div className="physics-box">
        <h4><Thermometer size={18} /> Physics of Temperature Distribution</h4>
        
        <div className="physics-highlight magenta">
          <strong>Governing Equation:</strong>
          <div className="equation-inline" style={{ display: 'block', marginTop: '0.5rem' }}>
            A₃·θ'' + A₁·Pr·Ec·(W')² + A₂·Pr·Ec·Ha²·W² = 0
          </div>
        </div>
        
        <p><strong>Heat Generation Sources:</strong></p>
        <ul>
          <li><strong>A₁·Pr·Ec·(W')²</strong> — Viscous dissipation (friction heating)</li>
          <li><strong>A₂·Pr·Ec·Ha²·W²</strong> — Joule heating (magnetic field effect)</li>
        </ul>
        
        <div className="physics-highlight gold">
          <strong>Effect of Eckert Number (Ec):</strong><br/>
          Increasing Ec → More viscous dissipation → Higher temperatures in the fluid.
          This represents the conversion of kinetic energy to thermal energy.
        </div>
        
        <div className="physics-highlight emerald">
          <strong>Effect of Biot Number (Bi):</strong><br/>
          Increasing Bi → Enhanced convective cooling at the upper plate → 
          Lower temperatures near η = 1.
        </div>
        
        <p><strong>Boundary Conditions:</strong></p>
        <ul>
          <li><strong>η = 0 (Lower plate):</strong> θ = 1 (fixed wall temperature)</li>
          <li><strong>η = 1 (Upper plate):</strong> θ' + Bi·θ = 0 (convective cooling)</li>
        </ul>
      </div>
    </div>
  );

  const renderEntropy = () => (
    <div className="visualization-section animate-slide-up">
      <ResultsPanel />
      
      <div className="chart-section">
        <div className="chart-header">
          <span className="dot gold"></span>
          <h3>Entropy Generation Components</h3>
        </div>
        <div className="chart-wrapper">
          <ResponsiveContainer width="100%" height="100%">
            <AreaChart data={solution.chartData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
              <defs>
                <linearGradient id="heatGrad" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#ff006e" stopOpacity={0.6}/>
                  <stop offset="95%" stopColor="#ff006e" stopOpacity={0.1}/>
                </linearGradient>
                <linearGradient id="fluidGrad" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#00d4ff" stopOpacity={0.6}/>
                  <stop offset="95%" stopColor="#00d4ff" stopOpacity={0.1}/>
                </linearGradient>
                <linearGradient id="magGrad" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#ffd700" stopOpacity={0.6}/>
                  <stop offset="95%" stopColor="#ffd700" stopOpacity={0.1}/>
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
              <XAxis dataKey="eta" stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 12 }} />
              <YAxis stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 12 }} />
              <Tooltip content={<CustomTooltip />} />
              <Area type="monotone" dataKey="Ns_heat" stackId="1" stroke="#ff006e" fill="url(#heatGrad)" name="Heat Transfer" />
              <Area type="monotone" dataKey="Ns_fluid" stackId="1" stroke="#00d4ff" fill="url(#fluidGrad)" name="Fluid Friction" />
              <Area type="monotone" dataKey="Ns_magnetic" stackId="1" stroke="#ffd700" fill="url(#magGrad)" name="Magnetic Field" />
            </AreaChart>
          </ResponsiveContainer>
        </div>
      </div>
      
      <div className="chart-section">
        <div className="chart-header">
          <span className="dot emerald"></span>
          <h3>Bejan Number Be(η)</h3>
        </div>
        <div className="chart-wrapper">
          <ResponsiveContainer width="100%" height="100%">
            <AreaChart data={solution.chartData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
              <defs>
                <linearGradient id="bejanGrad" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#00ff9f" stopOpacity={0.5}/>
                  <stop offset="95%" stopColor="#00ff9f" stopOpacity={0.05}/>
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
              <XAxis dataKey="eta" stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 12 }} />
              <YAxis domain={[0, 1]} stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 12 }} />
              <Tooltip content={<CustomTooltip />} />
              <Area type="monotone" dataKey="Be" stroke="#00ff9f" strokeWidth={2} fill="url(#bejanGrad)" name="Bejan Number" />
              {/* Reference line at 0.5 */}
              <Line type="monotone" data={[{eta: 0, ref: 0.5}, {eta: 1, ref: 0.5}]} dataKey="ref" stroke="rgba(255,255,255,0.3)" strokeDasharray="5 5" dot={false} />
            </AreaChart>
          </ResponsiveContainer>
        </div>
      </div>
      
      <div className="physics-box">
        <h4><BarChart3 size={18} /> Understanding Entropy Generation</h4>
        
        <p>Entropy generation measures thermodynamic irreversibility - energy that cannot be recovered.</p>
        
        <div className="physics-highlight">
          <strong>Total Entropy:</strong>
          <div className="equation-inline" style={{ display: 'block', marginTop: '0.5rem' }}>
            Ns = Ns,heat + Ns,fluid + Ns,magnetic
          </div>
        </div>
        
        <p><strong>Three Sources of Irreversibility:</strong></p>
        <ul>
          <li><strong style={{color: '#ff006e'}}>Heat Transfer (Ns,heat):</strong> Due to temperature gradients</li>
          <li><strong style={{color: '#00d4ff'}}>Fluid Friction (Ns,fluid):</strong> Due to viscous shear</li>
          <li><strong style={{color: '#ffd700'}}>Magnetic Field (Ns,magnetic):</strong> Due to Joule heating</li>
        </ul>
        
        <div className="physics-highlight emerald">
          <strong>Bejan Number (Be):</strong>
          <div className="equation-inline" style={{ display: 'block', marginTop: '0.5rem' }}>
            Be = Ns,heat / Ns,total
          </div>
          <p style={{ marginTop: '0.5rem', marginBottom: 0 }}>
            • <strong>Be {'>'} 0.5:</strong> Heat transfer irreversibility dominates<br/>
            • <strong>Be {'<'} 0.5:</strong> Friction + magnetic irreversibility dominates
          </p>
        </div>
        
        <p style={{ marginTop: '1rem' }}>
          <strong>Current Average Bejan Number:</strong> 
          <span style={{ color: 'var(--accent-emerald)', fontFamily: 'JetBrains Mono', marginLeft: '0.5rem' }}>
            {solution.avgBe.toFixed(4)}
          </span>
          {solution.avgBe > 0.5 ? 
            <span style={{ color: 'var(--accent-magenta)' }}> (Heat transfer dominated)</span> : 
            <span style={{ color: 'var(--accent-cyan)' }}> (Friction/Magnetic dominated)</span>
          }
        </p>
      </div>
    </div>
  );

  const renderVideos = () => (
    <div className="gallery-grid animate-slide-up">
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
            <p>Upload your video here</p>
          </div>
          <div className="video-info">
            <h3>{video.title}</h3>
            <p>{video.desc}</p>
          </div>
        </div>
      ))}
    </div>
  );

  const renderFigures = () => (
    <div className="gallery-grid animate-slide-up">
      {[
        { title: "Grid Convergence Study", desc: "Spectral accuracy validation" },
        { title: "Analytical Validation", desc: "Comparison with exact solutions" },
        { title: "Hartmann Number Effects", desc: "Velocity profiles for varying Ha" },
        { title: "Reynolds Number Study", desc: "Flow characteristics vs Re" },
        { title: "Prandtl Number Analysis", desc: "Thermal response to Pr" },
        { title: "Eckert Number Effects", desc: "Viscous dissipation impact" },
        { title: "3D Ha-Re Surface", desc: "Parameter space visualization" },
        { title: "3D Pr-Ec Surface", desc: "Thermal parameter coupling" },
        { title: "Entropy Distribution", desc: "Irreversibility sources" },
        { title: "Bejan Number Profiles", desc: "Dominance analysis" },
        { title: "Nanofluid Enhancement", desc: "Comparison with base fluid" },
        { title: "Literature Validation", desc: "Kigodi et al. (2025) comparison" }
      ].map((figure, index) => (
        <div key={index} className="gallery-item">
          <div className="gallery-item-image">
            <div className="gallery-item-placeholder"><Image /></div>
          </div>
          <div className="gallery-item-info">
            <h3>Fig {index + 1}: {figure.title}</h3>
            <p>{figure.desc}</p>
          </div>
        </div>
      ))}
    </div>
  );

  const renderTheory = () => (
    <div className="theory-content animate-slide-up">
      <div className="equation-card">
        <h3><Activity size={20} /> Momentum Equation</h3>
        <div className="equation">A₁·W'' - A₂·Ha²·W + G = 0</div>
        <p className="equation-description">
          Describes velocity distribution with viscous effects (A₁), magnetic damping (Ha), and pressure gradient (G).
        </p>
      </div>
      
      <div className="equation-card">
        <h3><Thermometer size={20} /> Energy Equation</h3>
        <div className="equation">A₃·θ'' + A₁·Pr·Ec·(W')² + A₂·Pr·Ec·Ha²·W² = 0</div>
        <p className="equation-description">
          Includes thermal conduction (A₃), viscous dissipation, and Joule heating from the magnetic field.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><Layers size={20} /> Boundary Conditions</h3>
        <div className="equation">
          η=0: W=0, θ=1<br/>
          η=1: W-λW'=Re, θ'+Bi·θ=0
        </div>
        <p className="equation-description">
          Lower plate is stationary with fixed temperature. Upper plate has slip and convective cooling.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><Gauge size={20} /> Engineering Quantities</h3>
        <div className="equation">
          Cf = A₁·dW/dη<br/>
          Nu = -A₃·dθ/dη
        </div>
        <p className="equation-description">
          Skin friction coefficient and Nusselt number measure wall shear and heat transfer.
        </p>
      </div>
      
      <div className="equation-card full-width">
        <h3><Droplets size={20} /> Nanofluid Property Correlations</h3>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem' }}>
          <div className="equation">ρnf = (1-φ)ρf + φρs</div>
          <div className="equation">μnf = μf/(1-φ)^2.5</div>
          <div className="equation">(ρCp)nf = (1-φ)(ρCp)f + φ(ρCp)s</div>
          <div className="equation">knf/kf = Maxwell model</div>
        </div>
        <p className="equation-description" style={{ marginTop: '1rem' }}>
          These correlations model effective nanofluid properties based on nanoparticle volume fraction φ.
        </p>
      </div>
    </div>
  );

  return (
    <div className="app">
      <header className="header">
        <div className="header-content">
          <div className="logo-section">
            <div className="logo-icon"><Zap size={24} /></div>
            <div className="logo-text">
              <h1>MHD Nanofluid Flow</h1>
              <p>Couette Flow Simulation</p>
            </div>
          </div>
          
          <nav className="nav-tabs">
            {[
              { id: 'simulation', icon: Play, label: 'Simulation' },
              { id: 'velocity', icon: Wind, label: 'Velocity' },
              { id: 'temperature', icon: Thermometer, label: 'Temperature' },
              { id: 'entropy', icon: BarChart3, label: 'Entropy' },
              { id: 'videos', icon: Video, label: 'Videos' },
              { id: 'figures', icon: Image, label: 'Figures' },
              { id: 'theory', icon: BookOpen, label: 'Theory' }
            ].map(tab => (
              <button
                key={tab.id}
                className={`nav-tab ${activeTab === tab.id ? 'active' : ''}`}
                onClick={() => setActiveTab(tab.id)}
              >
                <tab.icon size={16} />
                <span className="tab-text">{tab.label}</span>
              </button>
            ))}
          </nav>
        </div>
      </header>
      
      <main className="main-content">
        {activeTab === 'simulation' && renderSimulation()}
        {activeTab === 'velocity' && renderVelocity()}
        {activeTab === 'temperature' && renderTemperature()}
        {activeTab === 'entropy' && renderEntropy()}
        {activeTab === 'videos' && renderVideos()}
        {activeTab === 'figures' && renderFigures()}
        {activeTab === 'theory' && renderTheory()}
      </main>
      
      <FloatingControls />
      
      <footer className="footer">
        <p>
          <strong>Research:</strong> Thermal & MHD Analysis of Nanofluid Couette Flow<br/>
          <strong>Candidate:</strong> Mr. S.I. Mosala | <strong>Supervisor:</strong> Prof. O.D. Makinde<br/>
          Nelson Mandela University | December 2025
        </p>
      </footer>
    </div>
  );
}

export default App;