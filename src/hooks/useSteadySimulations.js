import { useState, useCallback } from 'react';
import { supabase, getUserId } from '../supabaseClient';

export const useSteadySimulations = () => {
  const [isLoading, setIsLoading] = useState(false);
  const [trainingData, setTrainingData] = useState([]);
  const [leaderboardData, setLeaderboardData] = useState([]);

  // Calculate quality score based on residuals
  const calculateQualityScore = useCallback((results) => {
    const momentumResidual = results.residualMomentum?.avgResidual || 0;
    const energyResidual = results.residualEnergy?.avgResidual || 0;
    
    // Lower residuals = higher quality
    const residualScore = 1 / (1 + momentumResidual + energyResidual);
    
    // Penalize extreme values
    const rangeScore = (results.maxW < 100 && results.maxTheta < 10) ? 1 : 0.5;
    
    return (residualScore * rangeScore) * 100;
  }, []);

  // Save simulation for ML training
  const saveSimulationForML = useCallback(async (params, results, optimizationGoal = null) => {
    try {
      const qualityScore = calculateQualityScore(results);
      
      const { data, error } = await supabase
        .from('steady_simulations')
        .insert({
          user_id: getUserId(),
          params: {
            Ha: params.Ha,
            Re: params.Re,
            Pr: params.Pr,
            Ec: params.Ec,
            Bi: params.Bi,
            lambda: params.lambda,
            G: params.G,
            A1: params.A1,
            A2: params.A2,
            A3: params.A3
          },
          results: {
            Cf_lower: results.Cf_lower,
            Cf_upper: results.Cf_upper,
            Nu_lower: results.Nu_lower,
            Nu_upper: results.Nu_upper,
            avgNs: results.avgNs,
            maxW: results.maxW,
            minTheta: results.minTheta,
            maxTheta: results.maxTheta,
            avgBe: results.avgBe
          },
          quality_score: qualityScore,
          optimization_goal: optimizationGoal,
          metadata: {
            userAgent: navigator.userAgent,
            timestamp: Date.now()
          }
        })
        .select()
        .single();

      if (error) throw error;
      return data;
    } catch (error) {
      console.error('Error saving simulation:', error);
      return null;
    }
  }, [calculateQualityScore]);

  // Load training data for ML
  const loadTrainingData = useCallback(async (limit = 1000) => {
    setIsLoading(true);
    try {
      const { data, error } = await supabase
        .from('steady_simulations')
        .select('*')
        .order('created_at', { ascending: false })
        .limit(limit);

      if (error) throw error;

      setTrainingData(data || []);
      return data;
    } catch (error) {
      console.error('Error loading training data:', error);
      return [];
    } finally {
      setIsLoading(false);
    }
  }, []);

  // Submit to leaderboard
  const submitToLeaderboard = useCallback(async (goal, params, results, score, username = null) => {
    try {
      const { data, error } = await supabase
        .from('optimization_leaderboard')
        .insert({
          user_id: getUserId(),
          username: username || `Researcher ${getUserId().slice(-4)}`,
          optimization_goal: goal,
          params: {
            Ha: params.Ha,
            Re: params.Re,
            Pr: params.Pr,
            Ec: params.Ec,
            Bi: params.Bi,
            lambda: params.lambda,
            G: params.G
          },
          results: {
            Cf_lower: results.Cf_lower,
            Nu_lower: results.Nu_lower,
            avgNs: results.avgNs,
            maxW: results.maxW
          },
          score: score
        })
        .select()
        .single();

      if (error) throw error;
      return data;
    } catch (error) {
      console.error('Error submitting to leaderboard:', error);
      return null;
    }
  }, []);

  // Get leaderboard
  const getLeaderboard = useCallback(async (goal, limit = 10) => {
    try {
      const { data, error } = await supabase
        .from('optimization_leaderboard')
        .select('*')
        .eq('optimization_goal', goal)
        .order('score', { ascending: false })
        .limit(limit);

      if (error) throw error;
      setLeaderboardData(data || []);
      return data || [];
    } catch (error) {
      console.error('Error loading leaderboard:', error);
      return [];
    }
  }, []);

  // Get statistics about training data
  const getTrainingStats = useCallback(async () => {
    try {
      const { data, error } = await supabase
        .from('steady_simulations')
        .select('id');

      if (error) throw error;

      return {
        totalSimulations: data?.length || 0,
        lastUpdated: new Date().toISOString()
      };
    } catch (error) {
      console.error('Error loading stats:', error);
      return { totalSimulations: 0, lastUpdated: null };
    }
  }, []);

  return {
    isLoading,
    trainingData,
    leaderboardData,
    saveSimulationForML,
    loadTrainingData,
    submitToLeaderboard,
    getLeaderboard,
    getTrainingStats
  };
};