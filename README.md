# WSN-Routing: Correlation-Aware Routing Analysis

Implementation and comparison of Wireless Sensor Network routing algorithms based on [this research paper](https://www.sciencedirect.com/science/article/abs/pii/S1570870511002320).

## 🚀 Quick Start

```bash
pip install networkx matplotlib numpy
python static_corr_routing.py
```

## 📊 Results

### Algorithm Comparison
<img width="4469" height="1769" alt="Image" src="https://github.com/user-attachments/assets/7ed0857e-5d3f-4a4c-ab1f-30cb8323c3a9" />

```mermaid
graph TD
    A[Algorithm Comparison Results] --> B[Energy Efficiency]
    A --> C[Execution Time]
    A --> D[Routing Trees]
    
    B --> B1[Brute Force: 4.63e-04 J/s]
    B --> B2[MER: 4.63e-04 J/s]
    B --> B3[CAR: 4.63e-04 J/s]
    
    C --> C1[Brute Force: 0.0174s]
    C --> C2[MER: 0.0010s - 18.1x faster]
    C --> C3[CAR: 0.0012s - 14.5x faster]
    
    D --> D1[MER Tree: Optimal Energy Path]
    D --> D2[CAR Tree: Correlation-Aware Path]
    
    style B1 fill:#ff9999
    style B2 fill:#99ccff
    style B3 fill:#99ff99
    style C2 fill:#99ccff
    style C3 fill:#99ff99
```

**Key Findings:**
- **MER (Minimum Energy Routing)**: 18-156x faster than brute force
- **CAR (Correlation-Aware Routing)**: Optimal solutions with correlation modeling
- **Both algorithms achieve optimal energy efficiency** (0% gap vs brute force)

### Scalability Analysis
```mermaid
graph LR
    A[Scalability Analysis] --> B[Network Size]
    A --> C[Performance Trends]
    
    B --> B1[3 Sources: 1.6x speedup]
    B --> B2[4 Sources: 3.4x speedup]
    B --> B3[5 Sources: 21.2x speedup]
    B --> B4[6 Sources: 156.1x speedup]
    B --> B5[7 Sources: CAR 2.59% better]
    
    C --> C1[Polynomial Complexity]
    C --> C2[Optimal Solutions]
    C --> C3[Correlation Benefits]
    
    style B1 fill:#e1f5fe
    style B2 fill:#e1f5fe
    style B3 fill:#e1f5fe
    style B4 fill:#e1f5fe
    style B5 fill:#c8e6c9
```

**Performance:**
- Excellent scalability across network sizes (3-7 sources)
- CAR shows advantages in larger networks (2.59% better than MER at 7 sources)
- Polynomial time complexity vs exponential brute force

### Network Topology
```mermaid
graph TD
    D[Sink<br/>D] --> Y1[Source Y1]
    D --> Y2[Source Y2]
    Y1 --> Y3[Source Y3]
    Y2 --> Y4[Source Y4]
    Y3 --> Y5[Source Y5]
    Y4 --> Y5
    Y1 --> Y2
    Y2 --> Y3
    Y3 --> Y4
    Y4 --> Y1
    Y5 --> Y1
    
    style D fill:#ffcccb
    style Y1 fill:#87ceeb
    style Y2 fill:#87ceeb
    style Y3 fill:#87ceeb
    style Y4 fill:#87ceeb
    style Y5 fill:#87ceeb
```

## 🏗️ Algorithms

1. **Brute Force**: Optimal but exponential complexity (≤6 sources)
2. **MER**: Game-theoretic approach minimizing total energy
3. **CAR**: Advanced algorithm considering spatial data correlations

## 📈 Network Topology

- **Pentagon Configuration**: Sources arranged around central sink 
- **Fully Connected**: All sources can communicate with each other tho two node can hop on sink
- **Visualization**: Blue nodes = sources, salmon node = sink, green arrows = optimal routes

## 🔬 Technical Details

- **Physical Layer**: Path loss model with SNR/BER calculations
- **Data Aggregation**: Sequential loss model with correlation function `q(d) = exp(-d/corr_scale)`
- **Game Theory**: Best-response dynamics for convergence to Nash equilibrium

---

**Note**: This implementation demonstrates that both MER and CAR achieve optimal solutions while providing significant computational advantages over brute force enumeration.

