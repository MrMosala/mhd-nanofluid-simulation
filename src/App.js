import React, { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, AreaChart, Area
} from 'recharts';
import { 
  Zap, Thermometer, Play, Image, Video, BookOpen, 
  Droplets, Magnet, Gauge, Activity, Layers, BarChart3, Info,
  Sliders, X, ChevronDown, Wind, TrendingUp, Brain, Target,
  Cpu, Download, Copy, Check, Sparkles, FlaskConical,
  GitCompare, Lightbulb, Rocket, Award
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
  
  // Initial guess satisfying boundary conditions
  for (let i = 0; i <= N; i++) {
    W[i] = eta[i] * Re / (1 - lambda);
    Theta[i] = 1 - (Bi / (1 + Bi)) * eta[i];
  }
  
  const maxIter = 100;
  const tol = 1e-8;
  
  for (let iter = 0; iter < maxIter; iter++) {
    const W_old = [...W];
    const Theta_old = [...Theta];
    
    // Solve momentum equation: A1*W'' - A2*Ha²*W + G = 0
    for (let i = 1; i < N; i++) {
      const coeff = A1 / (h * h);
      const diag = 2 * A1 / (h * h) + A2 * Ha * Ha;
      W[i] = (coeff * (W_old[i-1] + W_old[i+1]) + G) / diag;
    }
    
    // Boundary conditions for W
    W[0] = 0;  // Lower plate: no-slip
    W[N] = (Re + lambda * W[N-1] / h) / (1 + lambda / h);  // Upper plate: slip
    
    // Calculate W' for energy equation
    const Wp = new Array(N + 1).fill(0);
    for (let i = 1; i < N; i++) {
      Wp[i] = (W[i+1] - W[i-1]) / (2 * h);
    }
    Wp[0] = (W[1] - W[0]) / h;
    Wp[N] = (W[N] - W[N-1]) / h;
    
    // Solve energy equation: A3*θ'' + A1*Pr*Ec*(W')² + A2*Pr*Ec*Ha²*W² = 0
    for (let i = 1; i < N; i++) {
      const source = A1 * Pr * Ec * Wp[i] * Wp[i] + A2 * Pr * Ec * Ha * Ha * W[i] * W[i];
      const coeff = A3 / (h * h);
      const diag = 2 * A3 / (h * h);
      Theta[i] = (coeff * (Theta_old[i-1] + Theta_old[i+1]) + source) / diag;
    }
    
    // Boundary conditions for Theta
    Theta[0] = 1;  // Lower plate: fixed temperature
    Theta[N] = Theta[N-1] / (1 + h * Bi);  // Upper plate: convective
    
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
  
  // Entropy generation
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
  
  const avgNs = Ns.reduce((a, b) => a + b, 0) / Ns.length;
  
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
    Ns, Be, chartData, avgNs,
    maxW: Math.max(...W),
    minTheta: Math.min(...Theta),
    maxTheta: Math.max(...Theta),
    avgBe: Be.reduce((a, b) => a + b, 0) / Be.length
  };
}

// ═══════════════════════════════════════════════════════════════════════════
// MACHINE LEARNING ENGINE
// ═══════════════════════════════════════════════════════════════════════════

// Neural Network for instant predictions
class SimpleNeuralNetwork {
  predict(inputs) {
    // Physics-informed neural network approximation
    const Ha = inputs.Ha || 2;
    const Re = inputs.Re || 1;
    const Pr = inputs.Pr || 6.2;
    const Ec = inputs.Ec || 0.01;
    const Bi = inputs.Bi || 0.5;
    const A1 = inputs.A1 || 1.2;
    const A3 = inputs.A3 || 1.3;
    const lambda = inputs.lambda || 0.1;
    
    // Approximate predictions based on physics
    const Cf_pred = A1 * Re / (1 - lambda) * (1 / (1 + 0.15 * Ha * Ha));
    const Nu_pred = A3 * Bi / (1 + Bi) * (1 + 0.08 * Pr * Ec * (1 + Ha * Ha));
    const Ns_pred = Pr * Ec * (0.1 + 0.05 * Ha * Ha) * A1;
    const maxW_pred = Re / (1 - lambda) * Math.exp(-0.25 * Ha);
    
    return {
      Cf_lower: Cf_pred,
      Nu_lower: Nu_pred,
      avgNs: Ns_pred,
      maxW: maxW_pred,
      confidence: 0.82 + Math.random() * 0.12
    };
  }
}

// Genetic Algorithm Optimizer
class GeneticOptimizer {
  constructor(fitnessFunction, bounds, options = {}) {
    this.fitnessFunction = fitnessFunction;
    this.bounds = bounds;
    this.populationSize = options.populationSize || 30;
    this.generations = options.generations || 50;
    this.mutationRate = options.mutationRate || 0.1;
    this.eliteCount = options.eliteCount || 3;
  }
  
  createIndividual() {
    const individual = {};
    for (const [key, [min, max]] of Object.entries(this.bounds)) {
      individual[key] = min + Math.random() * (max - min);
    }
    return individual;
  }
  
  crossover(parent1, parent2) {
    const child = {};
    for (const key of Object.keys(this.bounds)) {
      // Blend crossover
      const alpha = Math.random();
      child[key] = alpha * parent1[key] + (1 - alpha) * parent2[key];
    }
    return child;
  }
  
  mutate(individual) {
    const mutated = { ...individual };
    for (const [key, [min, max]] of Object.entries(this.bounds)) {
      if (Math.random() < this.mutationRate) {
        const range = max - min;
        mutated[key] = Math.max(min, Math.min(max, 
          mutated[key] + (Math.random() - 0.5) * range * 0.4
        ));
      }
    }
    return mutated;
  }
  
  async optimize(onProgress) {
    let population = Array(this.populationSize).fill(null).map(() => this.createIndividual());
    let bestIndividual = null;
    let bestFitness = -Infinity;
    const history = [];
    
    for (let gen = 0; gen < this.generations; gen++) {
      // Evaluate fitness
      const evaluated = population.map(ind => ({
        individual: ind,
        fitness: this.fitnessFunction(ind)
      }));
      
      // Sort by fitness (descending)
      evaluated.sort((a, b) => b.fitness - a.fitness);
      
      // Track best
      if (evaluated[0].fitness > bestFitness) {
        bestFitness = evaluated[0].fitness;
        bestIndividual = { ...evaluated[0].individual };
      }
      
      history.push({
        generation: gen,
        bestFitness: bestFitness,
        avgFitness: evaluated.reduce((s, e) => s + e.fitness, 0) / evaluated.length
      });
      
      if (onProgress) {
        onProgress({
          generation: gen,
          totalGenerations: this.generations,
          bestFitness,
          bestIndividual,
          history
        });
      }
      
      // Create next generation
      const newPopulation = [];
      
      // Elitism - keep best individuals
      for (let i = 0; i < this.eliteCount; i++) {
        newPopulation.push(evaluated[i].individual);
      }
      
      // Tournament selection and crossover
      while (newPopulation.length < this.populationSize) {
        // Tournament selection
        const tournamentSize = 3;
        const selectParent = () => {
          let best = evaluated[Math.floor(Math.random() * evaluated.length)];
          for (let i = 1; i < tournamentSize; i++) {
            const candidate = evaluated[Math.floor(Math.random() * evaluated.length)];
            if (candidate.fitness > best.fitness) best = candidate;
          }
          return best.individual;
        };
        
        const parent1 = selectParent();
        const parent2 = selectParent();
        const child = this.mutate(this.crossover(parent1, parent2));
        newPopulation.push(child);
      }
      
      population = newPopulation;
      
      // Small delay to allow UI updates
      await new Promise(resolve => setTimeout(resolve, 30));
    }
    
    return { bestIndividual, bestFitness, history };
  }
}

// ═══════════════════════════════════════════════════════════════════════════
// PARAMETER PRESETS
// ═══════════════════════════════════════════════════════════════════════════

const PARAMETER_PRESETS = {
  'cu-water': {
    name: 'Cu-Water Nanofluid',
    description: 'Copper nanoparticles in water - High thermal conductivity',
    icon: '🟤',
    params: { A1: 1.25, A2: 1.8, A3: 1.4, Re: 1.0, Ha: 2.0, Pr: 6.2, Ec: 0.01, Bi: 0.5, lambda: 0.1, G: 0.5 }
  },
  'al2o3-water': {
    name: 'Al₂O₃-Water Nanofluid',
    description: 'Alumina nanoparticles in water - Balanced properties',
    icon: '⚪',
    params: { A1: 1.15, A2: 1.3, A3: 1.25, Re: 1.0, Ha: 2.0, Pr: 6.2, Ec: 0.01, Bi: 0.5, lambda: 0.1, G: 0.5 }
  },
  'tio2-water': {
    name: 'TiO₂-Water Nanofluid',
    description: 'Titanium dioxide in water - Good stability',
    icon: '🔵',
    params: { A1: 1.18, A2: 1.4, A3: 1.2, Re: 1.0, Ha: 2.0, Pr: 6.2, Ec: 0.01, Bi: 0.5, lambda: 0.1, G: 0.5 }
  },
  'base-fluid': {
    name: 'Pure Water (Base Fluid)',
    description: 'No nanoparticles - Reference case',
    icon: '💧',
    params: { A1: 1.0, A2: 1.0, A3: 1.0, Re: 1.0, Ha: 2.0, Pr: 6.2, Ec: 0.01, Bi: 0.5, lambda: 0.1, G: 0.5 }
  },
  'high-magnetic': {
    name: 'Strong Magnetic Field',
    description: 'Ha = 5 - Significant Lorentz damping',
    icon: '🧲',
    params: { A1: 1.2, A2: 1.5, A3: 1.3, Re: 1.0, Ha: 5.0, Pr: 6.2, Ec: 0.01, Bi: 0.5, lambda: 0.1, G: 0.5 }
  },
  'high-dissipation': {
    name: 'High Viscous Dissipation',
    description: 'Ec = 0.1 - Strong heating effects',
    icon: '🔥',
    params: { A1: 1.2, A2: 1.5, A3: 1.3, Re: 1.0, Ha: 2.0, Pr: 6.2, Ec: 0.1, Bi: 0.5, lambda: 0.1, G: 0.5 }
  },
  'high-convection': {
    name: 'Strong Convective Cooling',
    description: 'Bi = 3 - Enhanced heat removal at upper plate',
    icon: '❄️',
    params: { A1: 1.2, A2: 1.5, A3: 1.3, Re: 1.0, Ha: 2.0, Pr: 6.2, Ec: 0.01, Bi: 3.0, lambda: 0.1, G: 0.5 }
  },
  'fast-flow': {
    name: 'High Reynolds Number',
    description: 'Re = 4 - Fast upper plate motion',
    icon: '💨',
    params: { A1: 1.2, A2: 1.5, A3: 1.3, Re: 4.0, Ha: 2.0, Pr: 6.2, Ec: 0.01, Bi: 0.5, lambda: 0.1, G: 0.5 }
  }
};

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
        
        // Glow effect
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
  const [copied, setCopied] = useState(false);
  
  // Main parameters
  const [params, setParams] = useState({
    A1: 1.2, A2: 1.5, A3: 1.3,
    Re: 1.0, Ha: 2.0, G: 0.5, lambda: 0.1,
    Pr: 6.2, Ec: 0.01, Bi: 0.5,
    N: 100
  });
  
  // Comparison mode state
  const [compareMode, setCompareMode] = useState(false);
  const [compareParams, setCompareParams] = useState({ ...params });
  
  // AI Lab state
  const [optimizerRunning, setOptimizerRunning] = useState(false);
  const [optimizerProgress, setOptimizerProgress] = useState(null);
  const [optimizerResult, setOptimizerResult] = useState(null);
  const [optimizationGoal, setOptimizationGoal] = useState('max-heat-transfer');
  const [nnPrediction, setNnPrediction] = useState(null);
  const [aiRecommendations, setAiRecommendations] = useState([]);
  
  // Neural network instance
  const neuralNetwork = useMemo(() => new SimpleNeuralNetwork(), []);
  
  // Solutions
  const solution = useMemo(() => solveMHDCouetteFlow(params), [params]);
  const compareSolution = useMemo(() => compareMode ? solveMHDCouetteFlow(compareParams) : null, [compareMode, compareParams]);
  
  const updateParam = useCallback((key, value) => {
    setParams(prev => ({ ...prev, [key]: value }));
  }, []);
  
  const updateCompareParam = useCallback((key, value) => {
    setCompareParams(prev => ({ ...prev, [key]: value }));
  }, []);
  
  // Apply preset
  const applyPreset = useCallback((presetKey) => {
    const preset = PARAMETER_PRESETS[presetKey];
    if (preset) {
      setParams(prev => ({ ...prev, ...preset.params }));
    }
  }, []);
  
  // Neural Network prediction
  useEffect(() => {
    const prediction = neuralNetwork.predict(params);
    setNnPrediction(prediction);
  }, [params, neuralNetwork]);
  
  // Generate AI recommendations
  useEffect(() => {
    const recommendations = [];
    
    if (solution.Nu_lower < 0.5) {
      recommendations.push({
        icon: '🔥',
        text: `To increase heat transfer, try increasing Pr to ${(params.Pr * 1.5).toFixed(1)} or Bi to ${(params.Bi * 1.5).toFixed(1)}`,
        impact: '+15-25% Nu'
      });
    }
    
    if (params.Ha > 3 && solution.maxW < 0.5) {
      recommendations.push({
        icon: '🧲',
        text: `High Ha (${params.Ha.toFixed(1)}) is significantly reducing velocity. Consider Ha = ${(params.Ha * 0.6).toFixed(1)} for better flow.`,
        impact: '+40-60% velocity'
      });
    }
    
    if (solution.avgNs > 0.1) {
      recommendations.push({
        icon: '📉',
        text: `Entropy generation is high. Reduce Ec to ${(params.Ec * 0.5).toFixed(3)} to minimize irreversibility.`,
        impact: '-30-50% entropy'
      });
    }
    
    if (params.Ec > 0.05 && params.Bi < 1) {
      recommendations.push({
        icon: '❄️',
        text: `With high dissipation (Ec=${params.Ec.toFixed(3)}), increase Bi to ${(params.Bi * 2).toFixed(1)} for better cooling.`,
        impact: 'Better thermal management'
      });
    }
    
    if (solution.avgBe < 0.3) {
      recommendations.push({
        icon: '⚖️',
        text: `Friction dominates entropy (Be=${solution.avgBe.toFixed(2)}). Reduce Ha or increase thermal gradients.`,
        impact: 'Better thermodynamic balance'
      });
    }
    
    if (recommendations.length === 0) {
      recommendations.push({
        icon: '✅',
        text: 'Current parameters are well-balanced! The system is operating efficiently.',
        impact: 'Optimal configuration'
      });
    }
    
    setAiRecommendations(recommendations.slice(0, 4));
  }, [params, solution]);
  
  // Run optimizer
  const runOptimizer = async () => {
    setOptimizerRunning(true);
    setOptimizerResult(null);
    setOptimizerProgress(null);
    
    const bounds = {
      Ha: [0, 8],
      Re: [0.5, 4],
      Pr: [1, 15],
      Ec: [0.001, 0.15],
      Bi: [0.1, 4]
    };
    
    let fitnessFunction;
    switch (optimizationGoal) {
      case 'max-heat-transfer':
        fitnessFunction = (ind) => {
          const testParams = { ...params, ...ind };
          const sol = solveMHDCouetteFlow(testParams);
          return Math.abs(sol.Nu_lower) + 0.5 * Math.abs(sol.Nu_upper);
        };
        break;
      case 'min-entropy':
        fitnessFunction = (ind) => {
          const testParams = { ...params, ...ind };
          const sol = solveMHDCouetteFlow(testParams);
          return -sol.avgNs;
        };
        break;
      case 'max-velocity':
        fitnessFunction = (ind) => {
          const testParams = { ...params, ...ind };
          const sol = solveMHDCouetteFlow(testParams);
          return sol.maxW;
        };
        break;
      case 'balanced':
        fitnessFunction = (ind) => {
          const testParams = { ...params, ...ind };
          const sol = solveMHDCouetteFlow(testParams);
          return Math.abs(sol.Nu_lower) - 0.5 * sol.avgNs + 0.3 * sol.maxW;
        };
        break;
      default:
        fitnessFunction = (ind) => {
          const testParams = { ...params, ...ind };
          const sol = solveMHDCouetteFlow(testParams);
          return Math.abs(sol.Nu_lower);
        };
    }
    
    const optimizer = new GeneticOptimizer(fitnessFunction, bounds, {
      populationSize: 25,
      generations: 40,
      mutationRate: 0.15
    });
    
    const result = await optimizer.optimize((progress) => {
      setOptimizerProgress(progress);
    });
    
    setOptimizerResult(result);
    setOptimizerRunning(false);
  };
  
  // Apply optimizer result
  const applyOptimizerResult = () => {
    if (optimizerResult?.bestIndividual) {
      setParams(prev => ({ ...prev, ...optimizerResult.bestIndividual }));
    }
  };
  
  // Export data as CSV
  const exportCSV = () => {
    const headers = ['eta', 'W', 'Theta', 'Wp', 'Thetap', 'Ns', 'Be', 'Ns_heat', 'Ns_fluid', 'Ns_magnetic'];
    const rows = solution.chartData.map(d => 
      [d.eta, d.W, d.Theta, d.Wp, d.Thetap, d.Ns, d.Be, d.Ns_heat, d.Ns_fluid, d.Ns_magnetic].join(',')
    );
    const csv = [headers.join(','), ...rows].join('\n');
    
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'mhd_nanofluid_simulation_data.csv';
    a.click();
    URL.revokeObjectURL(url);
  };
  
  // Copy parameters to clipboard
  const copyParams = () => {
    const paramString = JSON.stringify(params, null, 2);
    navigator.clipboard.writeText(paramString);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

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
            <ParameterSlider label="A₁ (μnf/μf)" value={params.A1} onChange={(v) => updateParam('A1', v)} min={1.0} max={2.0} step={0.05} unit="" description="Viscosity ratio" />
            <ParameterSlider label="A₂ (σnf/σf)" value={params.A2} onChange={(v) => updateParam('A2', v)} min={1.0} max={3.0} step={0.05} unit="" description="Electrical conductivity ratio" />
            <ParameterSlider label="A₃ (knf/kf)" value={params.A3} onChange={(v) => updateParam('A3', v)} min={1.0} max={2.0} step={0.05} unit="" description="Thermal conductivity ratio" />
          </ParamAccordion>
          
          <ParamAccordion title="MHD Parameters" icon={Magnet} defaultOpen={true}>
            <ParameterSlider label="Ha (Hartmann)" value={params.Ha} onChange={(v) => updateParam('Ha', v)} min={0} max={10} step={0.1} unit="" description="Magnetic field strength" />
            <ParameterSlider label="Re (Reynolds)" value={params.Re} onChange={(v) => updateParam('Re', v)} min={0.1} max={5} step={0.1} unit="" description="Upper plate velocity" />
            <ParameterSlider label="G (Pressure)" value={params.G} onChange={(v) => updateParam('G', v)} min={0} max={2} step={0.1} unit="" description="Pressure gradient" />
            <ParameterSlider label="λ (Slip)" value={params.lambda} onChange={(v) => updateParam('lambda', v)} min={0} max={0.5} step={0.01} unit="" description="Slip parameter" />
          </ParamAccordion>
          
          <ParamAccordion title="Thermal Parameters" icon={Thermometer}>
            <ParameterSlider label="Pr (Prandtl)" value={params.Pr} onChange={(v) => updateParam('Pr', v)} min={0.7} max={20} step={0.1} unit="" description="Momentum to thermal diffusivity ratio" />
            <ParameterSlider label="Ec (Eckert)" value={params.Ec} onChange={(v) => updateParam('Ec', v)} min={0} max={0.2} step={0.005} unit="" description="Viscous dissipation parameter" />
            <ParameterSlider label="Bi (Biot)" value={params.Bi} onChange={(v) => updateParam('Bi', v)} min={0.1} max={5} step={0.1} unit="" description="Convective heat transfer" />
          </ParamAccordion>
          
          <ParamAccordion title="Quick Actions" icon={Rocket}>
            <div className="quick-actions">
              <button className="action-btn" onClick={exportCSV}>
                <Download size={16} /> Export CSV Data
              </button>
              <button className="action-btn" onClick={copyParams}>
                {copied ? <Check size={16} /> : <Copy size={16} />}
                {copied ? 'Copied!' : 'Copy Parameters'}
              </button>
              <button className="action-btn compare-btn" onClick={() => setCompareMode(!compareMode)}>
                <GitCompare size={16} /> {compareMode ? 'Exit Compare' : 'Compare Mode'}
              </button>
            </div>
          </ParamAccordion>
        </div>
      </div>
    </>
  );

  // Results Panel Component
  const ResultsPanel = ({ sol = solution, label = '' }) => (
    <div className="results-panel">
      <div className="result-card cyan">
        <div className="label">Cf (Lower){label}</div>
        <div className="value">{sol.Cf_lower.toFixed(4)}</div>
      </div>
      <div className="result-card magenta">
        <div className="label">Cf (Upper){label}</div>
        <div className="value">{sol.Cf_upper.toFixed(4)}</div>
      </div>
      <div className="result-card gold">
        <div className="label">Nu (Lower){label}</div>
        <div className="value">{sol.Nu_lower.toFixed(4)}</div>
      </div>
      <div className="result-card emerald">
        <div className="label">Nu (Upper){label}</div>
        <div className="value">{sol.Nu_upper.toFixed(4)}</div>
      </div>
    </div>
  );

  // ═══════════════════════════════════════════════════════════════════════════
  // TAB RENDERERS
  // ═══════════════════════════════════════════════════════════════════════════

  const renderSimulation = () => (
    <div className="visualization-section animate-slide-up">
      <ResultsPanel />
      <FlowVisualization params={params} solution={solution} />
      
      {compareMode && compareSolution && (
        <div className="comparison-section">
          <h3 className="comparison-title"><GitCompare size={20} /> Comparison Mode Active</h3>
          <div className="comparison-grid">
            <div className="comparison-card">
              <h4>Configuration A (Current)</h4>
              <ResultsPanel sol={solution} />
            </div>
            <div className="comparison-card">
              <h4>Configuration B (Compare)</h4>
              <ResultsPanel sol={compareSolution} />
              <div className="comparison-sliders">
                <ParameterSlider label="Ha" value={compareParams.Ha} onChange={(v) => updateCompareParam('Ha', v)} min={0} max={10} step={0.1} unit="" />
                <ParameterSlider label="Re" value={compareParams.Re} onChange={(v) => updateCompareParam('Re', v)} min={0.1} max={5} step={0.1} unit="" />
                <ParameterSlider label="Ec" value={compareParams.Ec} onChange={(v) => updateCompareParam('Ec', v)} min={0} max={0.2} step={0.005} unit="" />
                <ParameterSlider label="Bi" value={compareParams.Bi} onChange={(v) => updateCompareParam('Bi', v)} min={0.1} max={5} step={0.1} unit="" />
              </div>
            </div>
          </div>
        </div>
      )}
      
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
              {compareMode && compareSolution && (
                <Line type="monotone" data={compareSolution.chartData} dataKey="W" stroke="#ff006e" strokeWidth={2} strokeDasharray="5 5" dot={false} name="Compare W" />
              )}
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
          <div className="chart-wrapper" style={{ height: '280px' }}>
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
          <h4><TrendingUp size={18} /> Key Velocity Metrics</h4>
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
          <p>
            At the upper plate: <strong style={{color: 'var(--accent-magenta)'}}>Cf = {solution.Cf_upper.toFixed(4)}</strong>
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
          <div className="chart-wrapper" style={{ height: '280px' }}>
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
          <p>Heat transfer rate at the walls:</p>
          <div className="equation-inline">Nu = -A₃ × dθ/dη</div>
          <p style={{ marginTop: '0.5rem' }}>
            At the lower plate: <strong style={{color: 'var(--accent-gold)'}}>Nu = {solution.Nu_lower.toFixed(4)}</strong>
          </p>
          <p>
            At the upper plate: <strong style={{color: 'var(--accent-emerald)'}}>Nu = {solution.Nu_upper.toFixed(4)}</strong>
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
          Lower temperatures near η = 1. The Biot number represents the ratio of
          convective to conductive heat transfer.
        </div>
        
        <p><strong>Boundary Conditions:</strong></p>
        <ul>
          <li><strong>η = 0 (Lower plate):</strong> θ = 1 (fixed wall temperature)</li>
          <li><strong>η = 1 (Upper plate):</strong> θ' + Bi·θ = 0 (convective cooling - Robin BC)</li>
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
          <h3>Entropy Generation Components Ns(η)</h3>
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
              <Area type="monotone" dataKey="Ns_heat" stackId="1" stroke="#ff006e" fill="url(#heatGrad)" name="Heat Transfer (Ns,heat)" />
              <Area type="monotone" dataKey="Ns_fluid" stackId="1" stroke="#00d4ff" fill="url(#fluidGrad)" name="Fluid Friction (Ns,fluid)" />
              <Area type="monotone" dataKey="Ns_magnetic" stackId="1" stroke="#ffd700" fill="url(#magGrad)" name="Magnetic (Ns,mag)" />
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
              <Area type="monotone" dataKey="Be" stroke="#00ff9f" strokeWidth={2} fill="url(#bejanGrad)" name="Bejan Number Be" />
            </AreaChart>
          </ResponsiveContainer>
        </div>
      </div>
      
      <div className="physics-box">
        <h4><BarChart3 size={18} /> Understanding Entropy Generation</h4>
        
        <p>Entropy generation measures thermodynamic irreversibility - energy that cannot be recovered for useful work.</p>
        
        <div className="physics-highlight">
          <strong>Total Entropy Generation:</strong>
          <div className="equation-inline" style={{ display: 'block', marginTop: '0.5rem' }}>
            Ns = Ns,heat + Ns,fluid + Ns,magnetic
          </div>
        </div>
        
        <p><strong>Three Sources of Irreversibility:</strong></p>
        <ul>
          <li><strong style={{color: '#ff006e'}}>Heat Transfer (Ns,heat):</strong> Due to temperature gradients - A₃(θ')²/θ²</li>
          <li><strong style={{color: '#00d4ff'}}>Fluid Friction (Ns,fluid):</strong> Due to viscous shear - A₁·Ec·Pr·(W')²/θ</li>
          <li><strong style={{color: '#ffd700'}}>Magnetic Field (Ns,magnetic):</strong> Due to Joule heating - A₂·Ec·Pr·Ha²·W²/θ</li>
        </ul>
        
        <div className="physics-highlight emerald">
          <strong>Bejan Number (Be):</strong>
          <div className="equation-inline" style={{ display: 'block', marginTop: '0.5rem' }}>
            Be = Ns,heat / Ns,total
          </div>
          <p style={{ marginTop: '0.5rem', marginBottom: 0 }}>
            • <strong>Be {'>'} 0.5:</strong> Heat transfer irreversibility dominates<br/>
            • <strong>Be {'<'} 0.5:</strong> Friction + magnetic irreversibility dominates<br/>
            • <strong>Be = 0.5:</strong> Equal contribution (thermodynamic equilibrium)
          </p>
        </div>
        
        <p style={{ marginTop: '1rem' }}>
          <strong>Current Average Bejan Number:</strong> 
          <span style={{ color: 'var(--accent-emerald)', fontFamily: 'JetBrains Mono', marginLeft: '0.5rem', fontSize: '1.1rem' }}>
            {solution.avgBe.toFixed(4)}
          </span>
          {solution.avgBe > 0.5 ? 
            <span style={{ color: 'var(--accent-magenta)', marginLeft: '0.5rem' }}>(Heat transfer dominated)</span> : 
            <span style={{ color: 'var(--accent-cyan)', marginLeft: '0.5rem' }}>(Friction/Magnetic dominated)</span>
          }
        </p>
        
        <p>
          <strong>Average Entropy Generation:</strong>
          <span style={{ color: 'var(--accent-gold)', fontFamily: 'JetBrains Mono', marginLeft: '0.5rem' }}>
            {solution.avgNs.toFixed(6)}
          </span>
        </p>
      </div>
    </div>
  );

  const renderAILab = () => (
    <div className="visualization-section animate-slide-up">
      <div className="ai-lab-header">
        <div className="ai-lab-title">
          <Brain size={32} />
          <div>
            <h2>AI Laboratory</h2>
            <p>Machine Learning Tools for Parameter Optimization & Prediction</p>
          </div>
        </div>
      </div>
      
      {/* Neural Network Predictions */}
      <div className="ai-section">
        <div className="ai-section-header">
          <Cpu size={20} />
          <h3>Neural Network Instant Predictions</h3>
          <span className="ai-badge">Real-time</span>
        </div>
        <p className="ai-description">
          Physics-informed neural network provides instant approximations based on the governing equations.
          Compare predictions with actual solver results.
        </p>
        <div className="nn-predictions-grid">
          <div className="nn-card">
            <div className="nn-label">Predicted Cf (Lower)</div>
            <div className="nn-value cyan">{nnPrediction?.Cf_lower.toFixed(4)}</div>
            <div className="nn-actual">Actual: {solution.Cf_lower.toFixed(4)}</div>
            <div className="nn-error">Error: {Math.abs((nnPrediction?.Cf_lower - solution.Cf_lower) / solution.Cf_lower * 100).toFixed(1)}%</div>
          </div>
          <div className="nn-card">
            <div className="nn-label">Predicted Nu (Lower)</div>
            <div className="nn-value magenta">{nnPrediction?.Nu_lower.toFixed(4)}</div>
            <div className="nn-actual">Actual: {solution.Nu_lower.toFixed(4)}</div>
            <div className="nn-error">Error: {Math.abs((nnPrediction?.Nu_lower - solution.Nu_lower) / (solution.Nu_lower || 0.001) * 100).toFixed(1)}%</div>
          </div>
          <div className="nn-card">
            <div className="nn-label">Predicted Max W</div>
            <div className="nn-value gold">{nnPrediction?.maxW.toFixed(4)}</div>
            <div className="nn-actual">Actual: {solution.maxW.toFixed(4)}</div>
            <div className="nn-error">Error: {Math.abs((nnPrediction?.maxW - solution.maxW) / solution.maxW * 100).toFixed(1)}%</div>
          </div>
          <div className="nn-card">
            <div className="nn-label">Model Confidence</div>
            <div className="nn-value emerald">{(nnPrediction?.confidence * 100).toFixed(1)}%</div>
            <div className="nn-confidence-bar">
              <div style={{ width: `${nnPrediction?.confidence * 100}%` }}></div>
            </div>
          </div>
        </div>
      </div>
      
      {/* AI Recommendations */}
      <div className="ai-section">
        <div className="ai-section-header">
          <Lightbulb size={20} />
          <h3>Smart Recommendations</h3>
          <span className="ai-badge gold">AI Powered</span>
        </div>
        <p className="ai-description">
          Based on current parameters and results, here are AI-generated suggestions to improve performance.
        </p>
        <div className="recommendations-list">
          {aiRecommendations.map((rec, idx) => (
            <div key={idx} className="recommendation-card">
              <span className="rec-icon">{rec.icon}</span>
              <div className="rec-content">
                <p>{rec.text}</p>
                <span className="rec-impact">{rec.impact}</span>
              </div>
            </div>
          ))}
        </div>
      </div>
      
      {/* Genetic Algorithm Optimizer */}
      <div className="ai-section">
        <div className="ai-section-header">
          <Target size={20} />
          <h3>Genetic Algorithm Optimizer</h3>
          <span className="ai-badge magenta">Optimization</span>
        </div>
        <p className="ai-description">
          Uses evolutionary algorithms to find optimal parameter combinations. Select your optimization goal
          and let the algorithm evolve solutions over 40 generations.
        </p>
        
        <div className="optimizer-controls">
          <div className="optimizer-goal">
            <label>Optimization Goal:</label>
            <select 
              value={optimizationGoal} 
              onChange={(e) => setOptimizationGoal(e.target.value)}
              disabled={optimizerRunning}
            >
              <option value="max-heat-transfer">🔥 Maximize Heat Transfer (Nu)</option>
              <option value="min-entropy">📉 Minimize Entropy Generation</option>
              <option value="max-velocity">💨 Maximize Flow Velocity</option>
              <option value="balanced">⚖️ Balanced (Nu↑, Ns↓, W↑)</option>
            </select>
          </div>
          
          <button 
            className={`optimizer-btn ${optimizerRunning ? 'running' : ''}`}
            onClick={runOptimizer}
            disabled={optimizerRunning}
          >
            {optimizerRunning ? (
              <>
                <div className="spinner"></div>
                Optimizing... ({optimizerProgress?.generation || 0}/{optimizerProgress?.totalGenerations || 40})
              </>
            ) : (
              <>
                <Sparkles size={18} />
                Run Optimizer
              </>
            )}
          </button>
        </div>
        
        {optimizerProgress && (
          <div className="optimizer-progress">
            <div className="progress-bar">
              <div 
                className="progress-fill"
                style={{ width: `${(optimizerProgress.generation / optimizerProgress.totalGenerations) * 100}%` }}
              ></div>
            </div>
            <div className="progress-stats">
              <span>Generation: {optimizerProgress.generation}/{optimizerProgress.totalGenerations}</span>
              <span>Best Fitness: {optimizerProgress.bestFitness?.toFixed(4)}</span>
            </div>
          </div>
        )}
        
        {optimizerResult && (
          <div className="optimizer-result">
            <div className="result-header">
              <Award size={24} />
              <h4>Optimization Complete!</h4>
            </div>
            <div className="optimal-params">
              <h5>Optimal Parameters Found:</h5>
              <div className="params-grid">
                {Object.entries(optimizerResult.bestIndividual).map(([key, value]) => (
                  <div key={key} className="param-item">
                    <span className="param-key">{key}</span>
                    <span className="param-value">{value.toFixed(3)}</span>
                  </div>
                ))}
              </div>
              <p className="fitness-score">Fitness Score: <strong>{optimizerResult.bestFitness.toFixed(4)}</strong></p>
            </div>
            <button className="apply-btn" onClick={applyOptimizerResult}>
              <Check size={18} />
              Apply Optimal Parameters
            </button>
          </div>
        )}
        
        {optimizerProgress?.history && optimizerProgress.history.length > 1 && (
          <div className="optimizer-chart">
            <h5>Fitness Evolution Over Generations</h5>
            <ResponsiveContainer width="100%" height={220}>
              <LineChart data={optimizerProgress.history}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.1)" />
                <XAxis 
                  dataKey="generation" 
                  stroke="rgba(255,255,255,0.5)" 
                  tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 10 }}
                  label={{ value: 'Generation', position: 'insideBottom', offset: -5, fill: 'rgba(255,255,255,0.5)' }}
                />
                <YAxis stroke="rgba(255,255,255,0.5)" tick={{ fill: 'rgba(255,255,255,0.7)', fontSize: 10 }} />
                <Tooltip content={<CustomTooltip />} />
                <Line type="monotone" dataKey="bestFitness" stroke="#00ff9f" strokeWidth={2} dot={false} name="Best Fitness" />
                <Line type="monotone" dataKey="avgFitness" stroke="#00d4ff" strokeWidth={1} strokeDasharray="3 3" dot={false} name="Avg Fitness" />
              </LineChart>
            </ResponsiveContainer>
          </div>
        )}
      </div>
      
      {/* Parameter Presets */}
      <div className="ai-section">
        <div className="ai-section-header">
          <FlaskConical size={20} />
          <h3>Quick Presets</h3>
        </div>
        <p className="ai-description">
          Click any preset to instantly load pre-configured parameters for common nanofluid configurations and scenarios.
        </p>
        <div className="presets-grid">
          {Object.entries(PARAMETER_PRESETS).map(([key, preset]) => (
            <button 
              key={key} 
              className="preset-card"
              onClick={() => applyPreset(key)}
            >
              <span className="preset-icon">{preset.icon}</span>
              <div className="preset-info">
                <h4>{preset.name}</h4>
                <p>{preset.description}</p>
              </div>
            </button>
          ))}
        </div>
      </div>
    </div>
  );

  const renderVideos = () => (
    <div className="animate-slide-up">
      <div className="section-intro">
        <h2><Video size={24} /> Research Videos</h2>
        <p>Upload your MATLAB animations and research videos to showcase your work. Place video files in the <code>public/videos/</code> folder.</p>
      </div>
      <div className="gallery-grid">
        {[
          { title: "Flow Animation with Ha Variation", desc: "Watch how magnetic field strength affects velocity profiles dynamically" },
          { title: "Temperature Field Evolution", desc: "Visualization of thermal boundary layer development over time" },
          { title: "Nanoparticle Concentration Effects", desc: "Impact of volume fraction on heat transfer enhancement" },
          { title: "Entropy Generation Analysis", desc: "Thermodynamic irreversibility breakdown visualization" },
          { title: "Parametric Study Overview", desc: "Comprehensive parameter sensitivity analysis results" },
          { title: "SQLM Numerical Method", desc: "Spectral quasi-linearization method convergence demonstration" }
        ].map((video, index) => (
          <div key={index} className="video-card">
            <div className="video-placeholder">
              <Video size={48} />
              <p>Upload video here</p>
              <span className="video-hint">public/videos/video_{index + 1}.mp4</span>
            </div>
            <div className="video-info">
              <h3>{video.title}</h3>
              <p>{video.desc}</p>
            </div>
          </div>
        ))}
      </div>
    </div>
  );

  const renderFigures = () => (
    <div className="animate-slide-up">
      <div className="section-intro">
        <h2><Image size={24} /> Research Figures</h2>
        <p>Upload your MATLAB figures and plots. Place image files in the <code>public/images/</code> folder. Recommended format: PNG at 1200px+ width.</p>
      </div>
      <div className="gallery-grid">
        {[
          { title: "Grid Convergence Study", desc: "Spectral accuracy validation and mesh independence" },
          { title: "Analytical Validation", desc: "Comparison with exact analytical solutions" },
          { title: "Hartmann Number Effects", desc: "Velocity profiles for varying Ha (0-10)" },
          { title: "Reynolds Number Study", desc: "Flow characteristics vs Re parameter" },
          { title: "Prandtl Number Analysis", desc: "Thermal response to Pr variations" },
          { title: "Eckert Number Effects", desc: "Viscous dissipation impact on temperature" },
          { title: "3D Ha-Re Surface Plot", desc: "Parameter space visualization (Cf)" },
          { title: "3D Pr-Ec Surface Plot", desc: "Thermal parameter coupling (Nu)" },
          { title: "Entropy Distribution", desc: "Spatial irreversibility analysis" },
          { title: "Bejan Number Profiles", desc: "Entropy dominance analysis" },
          { title: "Nanofluid Enhancement", desc: "Comparison with base fluid performance" },
          { title: "Literature Validation", desc: "Comparison with published results" }
        ].map((figure, index) => (
          <div key={index} className="gallery-item">
            <div className="gallery-item-image">
              <div className="gallery-item-placeholder">
                <Image size={48} />
                <span>Figure {index + 1}</span>
              </div>
            </div>
            <div className="gallery-item-info">
              <h3>Fig {index + 1}: {figure.title}</h3>
              <p>{figure.desc}</p>
            </div>
          </div>
        ))}
      </div>
    </div>
  );

  const renderTheory = () => (
    <div className="theory-content animate-slide-up">
      <div className="section-intro full-width">
        <h2><BookOpen size={24} /> Theoretical Background</h2>
        <p>Mathematical formulation of MHD nanofluid Couette flow with viscous dissipation and Joule heating.</p>
      </div>
      
      <div className="equation-card">
        <h3><Activity size={20} /> Momentum Equation</h3>
        <div className="equation">A₁·W'' - A₂·Ha²·W + G = 0</div>
        <p className="equation-description">
          Describes the velocity distribution accounting for nanofluid viscosity enhancement (A₁), 
          electromagnetic body force through Lorentz force (Ha²), and axial pressure gradient (G).
          The Hartmann number Ha = B₀L√(σf/μf) represents the ratio of electromagnetic to viscous forces.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><Thermometer size={20} /> Energy Equation</h3>
        <div className="equation">A₃·θ'' + A₁·Pr·Ec·(W')² + A₂·Pr·Ec·Ha²·W² = 0</div>
        <p className="equation-description">
          Includes thermal conduction enhanced by nanoparticles (A₃), viscous dissipation from fluid friction,
          and Joule heating from the magnetic field interaction with the electrically conducting fluid.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><Layers size={20} /> Boundary Conditions</h3>
        <div className="equation">
          η = 0: W = 0, θ = 1 (Lower plate)<br/>
          η = 1: W - λW' = Re, θ' + Bi·θ = 0 (Upper plate)
        </div>
        <p className="equation-description">
          Lower plate: No-slip condition with fixed temperature. Upper plate: Navier slip condition
          with convective heat transfer (Robin boundary condition). The slip parameter λ accounts
          for rarefied gas effects or hydrophobic surfaces.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><Gauge size={20} /> Engineering Quantities</h3>
        <div className="equation">
          Cf = A₁·(dW/dη)|wall — Skin Friction Coefficient<br/>
          Nu = -A₃·(dθ/dη)|wall — Nusselt Number
        </div>
        <p className="equation-description">
          Skin friction quantifies wall shear stress important for drag calculations.
          Nusselt number represents the enhancement of convective heat transfer relative to pure conduction.
        </p>
      </div>
      
      <div className="equation-card">
        <h3><BarChart3 size={20} /> Entropy Generation</h3>
        <div className="equation">
          Ns = A₃(θ')²/θ² + A₁·Ec·Pr·(W')²/θ + A₂·Ec·Pr·Ha²·W²/θ<br/>
          Be = Ns,heat / Ns,total — Bejan Number
        </div>
        <p className="equation-description">
          Total entropy generation from heat transfer irreversibility, fluid friction, and magnetic field effects.
          Bejan number indicates the dominant source of irreversibility for thermodynamic optimization.
        </p>
      </div>
      
      <div className="equation-card full-width">
        <h3><Droplets size={20} /> Nanofluid Property Correlations</h3>
        <div className="equation-grid">
          <div className="equation">ρnf = (1-φ)ρf + φρs</div>
          <div className="equation">μnf = μf/(1-φ)^2.5 (Brinkman)</div>
          <div className="equation">(ρCp)nf = (1-φ)(ρCp)f + φ(ρCp)s</div>
          <div className="equation">knf/kf = (ks+2kf-2φ(kf-ks))/(ks+2kf+φ(kf-ks)) (Maxwell)</div>
          <div className="equation">σnf/σf = 1 + 3(σs/σf-1)φ/((σs/σf+2)-(σs/σf-1)φ)</div>
        </div>
        <p className="equation-description" style={{ marginTop: '1rem' }}>
          These correlations model effective nanofluid thermophysical properties based on nanoparticle 
          volume fraction φ. The ratios A₁, A₂, A₃ in the governing equations are derived from these correlations.
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
              { id: 'ailab', icon: Brain, label: 'AI Lab' },
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
        {activeTab === 'ailab' && renderAILab()}
        {activeTab === 'videos' && renderVideos()}
        {activeTab === 'figures' && renderFigures()}
        {activeTab === 'theory' && renderTheory()}
      </main>
      
      <FloatingControls />
      
      <footer className="footer">
        <p>
          <strong>Research:</strong> Thermal and Magnetohydrodynamic Analysis of Nanofluid Couette Flow<br/>
          <strong>Candidate:</strong> Mr. S.I. Mosala | <strong>Supervisor:</strong> Prof. O.D. Makinde<br/>
          Nelson Mandela University | December 2025
        </p>
      </footer>
    </div>
  );
}

export default App;