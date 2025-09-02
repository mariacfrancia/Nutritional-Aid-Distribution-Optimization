# Nutritional Aid Distribution Optimization

## Overview
Optimized **distribution of nutritional aid** to villages with acute malnutrition using **GAMS**. Focused on **facility location, transportation cost, equity, and compromise solutions**.  

**Skills & Keywords:** Optimization, MIP, Facility Location, Logistics, Resource Allocation, Multi-Criteria Decision-Making, Pareto Frontiers.

---

## Problem
Plan food distribution from **Massinga** to rural villages:  

- Minimize **number of centers** (max 4).  
- Minimize **total distance traveled** by population (weighted by food).  
- Minimize **distribution costs**: center operation, empty vehicle travel, cargo transport.  
- Ensure **equity** under 80% supply.  
- Provide **compromise solutions** balancing distance, cost, and equity.

---

## Models

1. **Minimize Centers:** Optimal locations (`X(i)`) and village assignments (`Y(i,j)`).  
2. **Minimize Distance:** Reduce travel distance (`Y(i,j)`).  
3. **Minimize Costs:** Optimal distribution plan (`Q(i,j)`, `V(i,j)`), center costs.  
4. **Equity under Partial Supply:** Fair allocation (`R(i,j)`), minimize deviation (`mdesv`).  
5. **Compromise Solutions:** Weighted combination of distance, cost, and equity for decision-making.

---

## Multi-Criteria Analysis
- **Payment Matrix:** Trade-offs between distance, cost, and equity.  
- **Pareto Frontiers:** 2D and 3D efficient solutions.  
- **Decision Support:** Stakeholder-friendly compromise solutions.

- GAMS Documentation: [https://www.gams.com](https://www.gams.com)  
- Facility Location & Multi-Criteria Optimization literature.
