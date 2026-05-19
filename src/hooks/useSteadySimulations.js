import { useState, useCallback } from 'react';
import { supabase, getUserId } from '../supabaseClient';

export const useSteadySimulations = () => {
  const [isLoading, setIsLoading] = useState(false);
  const [trainingData, setTrainingData] = useState([]);
  const [leaderboardData, setLeaderboardData] = useState([]);

  const saveSimulationForML = useCallback(async (params, results, optimizationGoal = null) => {
    try {
      const { data, error } = await supabase
        .from('steady_simulations')
        .insert({
          user_id: getUserId(),
          params: {
            Ha: params.Ha, Re: params.Re, Pr: params.Pr,
            Ec: params.Ec, Bi: params.Bi, lambda: params.lambda,
            G: params.G, A1: params.A1, A2: params.A2, A3: params.A3
          },
          results: {
            Cf_lower: results.Cf_lower, Cf_upper: results.Cf_upper,
            Nu_lower: results.Nu_lower, Nu_upper: results.Nu_upper,
            avgNs: results.avgNs, maxW: results.maxW,
            minTheta: results.minTheta, maxTheta: results.maxTheta,
            avgBe: results.avgBe
          },
          optimization_goal: optimizationGoal
        })
        .select()
        .single();

      if (error) throw error;
      return data;
    } catch (error) {
      console.error('Error saving simulation:', error);
      return null;
    }
  }, []);

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

  const submitToLeaderboard = useCallback(async (goal, params, results, score, username = null) => {
    try {
      const { data, error } = await supabase
        .from('leaderboard')
        .insert({
          user_id: getUserId(),
          username: username || `Researcher ${getUserId().slice(-4)}`,
          optimization_goal: goal,
          params: {
            Ha: params.Ha, Re: params.Re, Pr: params.Pr,
            Ec: params.Ec, Bi: params.Bi, lambda: params.lambda, G: params.G
          },
          results: {
            Cf_lower: results.Cf_lower, Nu_lower: results.Nu_lower,
            avgNs: results.avgNs, maxW: results.maxW
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

  const getLeaderboard = useCallback(async (goal, limit = 10) => {
    try {
      const { data, error } = await supabase
        .from('leaderboard')
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
    isLoading, trainingData, leaderboardData,
    saveSimulationForML, loadTrainingData,
    submitToLeaderboard, getLeaderboard, getTrainingStats
  };
};