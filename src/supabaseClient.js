import { createClient } from '@supabase/supabase-js';

const supabaseUrl = process.env.REACT_APP_SUPABASE_URL || 'https://nmnjcfuvymuvztsqqiyr.supabase.co';
const supabaseAnonKey = process.env.REACT_APP_SUPABASE_ANON_KEY || 'sb_publishable_eAuH0gpBz1umvvpO8Wf5mw_Lg91gLs3';

export const supabase = createClient(supabaseUrl, supabaseAnonKey, {
  realtime: {
    params: {
      eventsPerSecond: 10
    }
  }
});

export const getUserId = () => {
  let userId = localStorage.getItem('mhd_user_id');
  if (!userId) {
    userId = `user_${Math.random().toString(36).substr(2, 9)}`;
    localStorage.setItem('mhd_user_id', userId);
  }
  return userId;
};

export const getRandomColor = () => {
  const colors = [
    '#8b5cf6', '#ec4899', '#fbbf24', '#14b8a6',
    '#06b6d4', '#10b981', '#f97316', '#6366f1',
  ];
  return colors[Math.floor(Math.random() * colors.length)];
};