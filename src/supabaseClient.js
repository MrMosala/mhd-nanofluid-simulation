import { createClient } from '@supabase/supabase-js';

const supabaseUrl = process.env.REACT_APP_SUPABASE_URL || 'https://zflrbhiembxgbjiehycn.supabase.co';
const supabaseAnonKey = process.env.REACT_APP_SUPABASE_ANON_KEY || 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InpmbHJiaGllbWJ4Z2JqaWVoeWNuIiwicm9sZSI6ImFub24iLCJpYXQiOjE3NjY2NDAzMzUsImV4cCI6MjA4MjIxNjMzNX0.5h6kDnUXXAljbUK705XW4E4FZ4LpEzTF8viTEwPz4FI';

export const supabase = createClient(supabaseUrl, supabaseAnonKey, {
  realtime: {
    params: {
      eventsPerSecond: 10
    }
  }
});

// Helper to generate user ID
export const getUserId = () => {
  let userId = localStorage.getItem('mhd_user_id');
  if (!userId) {
    userId = `user_${Math.random().toString(36).substr(2, 9)}`;
    localStorage.setItem('mhd_user_id', userId);
  }
  return userId;
};

// Helper to generate random color
export const getRandomColor = () => {
  const colors = [
    '#8b5cf6', // purple
    '#ec4899', // pink
    '#fbbf24', // amber
    '#14b8a6', // teal
    '#06b6d4', // cyan
    '#10b981', // emerald
    '#f97316', // orange
    '#6366f1', // indigo
  ];
  return colors[Math.floor(Math.random() * colors.length)];
};