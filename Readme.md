# Parallel N-body simulations of planetary systems: a direct approach
https://arxiv.org/pdf/2208.13562.pdf


Direct gravitational simulations of n-body systems have a time complexity O(n2), which gets computation- ally expensive as the number of bodies increases. Distributing this workload to multiple cores significantly speeds up the computation and is the fundamental principle behind parallel computing. This project involves simulating (evolving) our solar system for the next 1000 years (from 2015 to 3015) on the BlueCrystal super- computer. The gravitational bodies (planets, moons, asteroids) were successfully simulated, and the initial states (mass, position and velocity vectors) of the solar system objects were obtained via NASA’s JPL Horizons web interface. Two parallel computing domains are investigated: shared and distributed memory systems.
