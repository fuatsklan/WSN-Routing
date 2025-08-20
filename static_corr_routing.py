"""
This code reproduce and may be improve following article:
https://www.sciencedirect.com/science/article/abs/pii/S1570870511002320

Some TO DO:
    - Real time simulation and implementation.
    - Larger network comparison.
"""

import networkx as nx
import matplotlib.pyplot as plt
import math
import heapq
from collections import defaultdict, deque
import itertools
import time

# Simple topology

def pentagon(
        sink_pos,
        pentagon_center=(5,5),
        radius=3,
        num_sources=5,
        sink_neighbors=2
    ):
    G = nx.DiGraph()
    pos = {}
    pos['D'] = sink_pos
    G.add_node('D')

    cx, cy = pentagon_center
    sources = []
    for i in range(num_sources):
        angle = -2 * math.pi * i / num_sources
        x = cx + radius * math.cos(angle)
        y = cy + radius * math.sin(angle)
        src = f'Y{i+1}'
        pos[src] = (x, y)
        G.add_node(src)
        sources.append(src)

    # fully connected among sources
    for i in sources:
        for j in sources:
            if i != j:
                G.add_edge(i, j)

    # connect first sink_neighbors sources to sink
    for i, src in enumerate(sources):
        if i < max(0, min(sink_neighbors, len(sources))):
            G.add_edge(src, 'D')

    return G, pos, sources


# PHY / Link metrics

def pc_to_gamma_target(M=1000, Pc_target = 0.999):
    inner = 1.0 - (Pc_target ** (1.0 / M))
    inner = max(1e-12, min(1-1e-12, inner))
    return -2.0 * math.log(inner)

def compute_link_metrics_tdma(G, pos,
                              P=0.05,   # 50 mW
                              W=2e5,    # 200 kHz
                              r2=50e-6, # noise power
                              alpha=2.2,# path-loss exponent
                              M=1000, m=800,
                              Pc_target=0.999):

    def distance(u, v):
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        return math.hypot(x1 - x2, y1 - y2)

    def channel_gain(d):
        return (d ** (-alpha)) if d > 0 else 1.0

    gamma_target = pc_to_gamma_target(M=M, Pc_target=Pc_target)

    for u, v in G.edges():
        d = distance(u, v)
        h = channel_gain(d)
        L_min = (gamma_target * r2) / (P * h) if (P * h) > 0 else float('inf')  # min time duration
        L = max(1.0, L_min)
        X = W / L                                                                # slot rate basically
        snr = (P * h / L) / r2 if r2 > 0 else float('inf')
        BER = 0.5 * math.exp(-0.5 * snr)
        Pc_bit = max(0.0, min(1.0, 1.0 - 2.0 * BER))
        Pc_pkt = (Pc_bit ** M) if Pc_bit > 0 else 0.0
        Eb = (M * P) / (m * X * Pc_pkt) if (X * Pc_pkt) > 1e-18 else float('inf')   # Expected energy per bit
        G.edges[u, v].update({
            'distance': d, 'h': h, 'L': L, 'X': X,
            'snr': snr, 'Pc_bit':Pc_bit, 'Pc_pkt': Pc_pkt,
            'Eb': Eb, 'BER': BER
        })

def print_edge_metrics(G):
    print("=== Edge Metrics ===")
    for u, v in G.edges():
        edge = G.edges[u, v]
        print(f"{u} → {v}")
        print(f"  Distance       : {edge.get('distance'):.2f} m")
        print(f"  Channel Gain h : {edge.get('h'):.5e}")
        print(f"  Path Loss L    : {edge.get('L'):.5e}")
        print(f"  SNR            : {edge.get('snr'):.5f}")
        print(f"  BER            : {edge.get('BER'):.5e}")
        print(f"  P_success (bit): {edge.get('Pc_bit'):.5f}")
        print(f"  P_success (pkt): {edge.get('Pc_pkt'):.5f}")
        print(f"  Rate X         : {edge.get('X'):.2f} bps")
        print(f"  Energy/bit Eb  : {edge.get('Eb'):.5e} J")
        print("")


# Agregation model

def compute_aggregated_rates(G, pos, parent, W_raw, alpha_f=0.5, corr_scale=3.0):
    """
    Sequential loss model
    Closest-first child merge order
    Forgetting factor applied if not raw data merge only
    Correlations decay with distance
    """
    def correlation_q(d, corr_scale=3.0):
        return max(0.0, min(1.0, math.exp(-d / float(corr_scale))))

    def aggregate_two(Wa, Wb, q_eff):
        return max(Wa, Wb) + (1.0 - q_eff) * min(Wa, Wb)

    # Build children map
    kids = defaultdict(list)
    for u, p in parent.items():
        if p is not None:
            kids[p].append(u)

    W_agg = {n:None for n in G.nodes()}
    is_agg =  {n: False for n in G.nodes()}

    def dfs(n):
        if n == 'D':
            W_agg[n] = 0.0
            is_agg[n] = True
            return W_agg[n], is_agg[n]

        if len(kids[n]) == 0:
            W_agg[n] = W_raw[n]
            is_agg[n] = False
            return W_agg[n], is_agg[n]

        # Internal nodes: sort children by higher correlation, then shorter distance
        cq_list = []
        for c in kids[n]:
            dx = pos[n][0] - pos[c][0]
            dy = pos[n][1] - pos[c][1]
            d = math.hypot(dx, dy)
            q_nc = correlation_q(d, corr_scale=corr_scale)
            cq_list.append((c, q_nc, d))
        cq_list.sort(key=lambda t: (-t[1], t[2], t[0]))

        Wa = W_raw[n]
        aggflag = False
        for (c, q_nc, _d) in cq_list:
            Wc, is_c_agg = dfs(c)
            q_eff = q_nc if (not aggflag and not is_c_agg) else (alpha_f * q_nc)
            Wa = aggregate_two(Wa, Wc, q_eff)
            aggflag = True

        W_agg[n] = Wa
        is_agg[n] = aggflag
        return W_agg[n], is_agg[n]

    for n in G.nodes():
        if W_agg[n] is None:
            dfs(n)
    return W_agg, is_agg


# Utilities

def path_to_sink(parent, start):
    path = []
    cur = start
    seen = {cur}
    while cur != 'D':
        p = parent.get(cur)
        if p is None or p in seen:
            return []
        path.append((cur, p))
        cur = p
        seen.add(cur)
    return path

def extract_routes(parent):
    routes = {}
    for n in parent:
        if n == 'D':
            continue
        path_nodes = [n]
        cur = n
        seen = {cur}
        while cur != 'D':
            p = parent.get(cur)
            if p is None or p in seen:
                path_nodes = []
                break
            path_nodes.append(p)
            cur = p
            seen.add(cur)
        if path_nodes:
            routes[n] = path_nodes
    return routes

def total_energy(G, parent, W_agg):
    tot = 0.0
    for u, p in parent.items():
        if p is None or u == 'D':
            continue
        Eb = G.edges[u, p]['Eb']
        if Eb == float('inf'):
            return float('inf')
        tot += Eb * W_agg[u]
    return tot

def potential_energy(G, parent, W_agg):
    """Exact potential for CAR: sum_f Eb(f->p)*W_agg[f]."""
    return total_energy(G, parent, W_agg)

def capacity_ok(G, parent, W_agg, B=None):
    if B is None:
        return True
    for u, p in parent.items():
        if p is None or u == 'D':
            continue
        X = G.edges[u, p]['X']   # capacity (bps)
        Wk = W_agg[u]            # offered traffic (bps)
        if Wk <= 0:
            continue
        rho = Wk / X
        if rho > B + 1e-12:
            return False
    return True

def plot_tree(G, pos, parent, title="Routing Tree"):
    plt.figure(figsize=(8, 8))
    sources = [n for n in G.nodes() if n != 'D']
    nx.draw_networkx_edges(G, pos, edgelist=G.edges(), alpha=0.15, edge_color='gray', arrowsize=12, arrowstyle='->')
    tree_edges = [(u, p) for u, p in parent.items() if p is not None]
    nx.draw_networkx_edges(G, pos, edgelist=tree_edges, width=2.6, edge_color='green', arrowsize=14, arrowstyle='->')
    nx.draw_networkx_nodes(G, pos, nodelist=sources, node_color='skyblue', node_size=500, label='Sources')
    nx.draw_networkx_nodes(G, pos, nodelist=['D'], node_color='salmon', node_size=700, label='Sink')
    nx.draw_networkx_labels(G, pos, font_size=12)
    plt.title(title)
    plt.axis('equal'); plt.axis('off')
    plt.legend()
    plt.show()

def plot_tree_on_axis(G, pos, parent, ax, title):
    sources = [n for n in G.nodes() if n != 'D']
    nx.draw_networkx_edges(G, pos, edgelist=G.edges(), alpha=0.15,
                           edge_color='gray', arrowsize=12, arrowstyle='->', ax=ax)
    tree_edges = [(u, p) for u, p in parent.items() if p is not None]
    nx.draw_networkx_edges(G, pos, edgelist=tree_edges, width=2.6,
                           edge_color='green', arrowsize=14, arrowstyle='->', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=sources, node_color='skyblue',
                           node_size=300, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=['D'], node_color='salmon',
                           node_size=400, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=10, ax=ax)
    ax.set_title(title, fontsize=12)
    ax.axis('equal')
    ax.axis('off')


# Brute Force solution

def enumerate_all_trees(G, pos, sources,
                        W_raw_val=200.0,
                        Pc_target=0.99,
                        alpha_f=0.5,
                        corr_scale=3.0,
                        B=None,
                        verbose=True):
    """
    Return a list sorted by total energy considering all possible valid parent vectors.
    """
    compute_link_metrics_tdma(G, pos, Pc_target=Pc_target)

    W_raw = {n: 0.0 for n in G.nodes()}
    for s in sources:
        W_raw[s] = W_raw_val
    W_raw['D'] = 0.0

    # candidate parent options
    parent_options = {}
    for u in sources:
        opts = []
        for v in G.successors(u):
            Eb = G.edges[u, v]['Eb']
            if Eb < float('inf'):
                opts.append(v)
        parent_options[u] = sorted(opts)

    if any(len(vs) == 0 for vs in parent_options.values()):
        if verbose:
            print("[enumerate] some nodes have no feasible outgoing edges.")
        return []

    results = []
    ordered_sources = sorted(sources)
    all_choices = itertools.product(*[parent_options[u] for u in ordered_sources])

    count_total = 0
    count_valid = 0
    for choice_tuple in all_choices:
        count_total += 1
        parent = {'D': None}
        for u, v in zip(ordered_sources, choice_tuple):
            parent[u] = v

        # check validity (reach sink, no cycles)
        valid = True
        for u in ordered_sources:
            path = path_to_sink(parent, u)
            if not path:
                valid = False
                break
        if not valid:
            continue

        # compute aggregated rates and total energy
        W_agg, _ = compute_aggregated_rates(G, pos, parent, W_raw, alpha_f=alpha_f,
                                            corr_scale=corr_scale)
        if not capacity_ok(G, parent, W_agg, B=B):
            continue
        total = total_energy(G, parent, W_agg)
        routes = extract_routes(parent)
        results.append({
            'energy': total,
            'parent': dict(parent),
            'routes': routes,
            'W_agg': W_agg,
        })
        count_valid += 1

    if verbose:
        print(f'[enumerate] Explored {count_total} parent assignments; valid trees: {count_valid}')

    results.sort(key = lambda d: d['energy'])
    return results


# Mer & Car Game Solvers

def build_routes_and_users(parent):
    routes = extract_routes(parent)
    users_at = {n: set() for n in parent.keys()}
    for s, path in routes.items():
        for node in path[:-1]:  # exclude D
            users_at[node].add(s)
    return routes, users_at

def facility_energy(G, parent, W_agg, f):
    p = parent.get(f)
    if p is None or f == 'D':
        return 0.0
    Eb = G.edges[f, p]['Eb']
    if Eb == float('inf'):
        return float('inf')
    return Eb * W_agg[f]

def path_to_root(parent, i):
    if i == 'D':
        return ['D']
    cur = i
    seen = set([cur])
    seq = [cur]
    while cur != 'D':
        p = parent.get(cur)
        if p is None or p in seen:
            return []
        seq.append(p)
        seen.add(p)
        cur = p
    return seq

def player_cost_MER(G, parent, routes, i):
    """u_i^MER = sum Eb along route(i)."""
    path = routes.get(i, [])
    if not path:
        return float('inf')
    tot = 0.0
    for f in path[:-1]:  # exclude D
        p = parent.get(f)
        if p is None:
            return float('inf')
        Eb = G.edges[f, p]['Eb']
        tot += Eb
    return tot

def player_cost_CAR_from_cached(G, parent, route_nodes, W_with, W_without):
    """Helper: compute CAR marginal energy on a given path using cached W_agg."""
    tot = 0.0
    for f in route_nodes[:-1]:
        wf_with = facility_energy(G, parent, W_with, f)
        wf_wo   = facility_energy(G, parent, W_without, f)
        if wf_with == float('inf'):
            return float('inf')
        delta = wf_with - wf_wo
        if delta > 0:
            tot += delta
    return tot

def creates_cycle(parent, u, v):
    cur = v
    seen = set()
    while cur is not None and cur not in seen:
        if cur == u:
            return True
        seen.add(cur)
        cur = parent.get(cur)
    return False

def dijkstra_mer_parent(G, sink='D'):
    """Compute MER tree (min-sum Eb to sink) via Dijkstra on reversed graph."""
    dist = {n: float('inf') for n in G.nodes()}
    dist[sink] = 0.0
    prev = {n: None for n in G.nodes()}  # predecessor toward sink

    pq = [(0.0, sink)]
    while pq:
        d, x = heapq.heappop(pq)
        if d != dist[x]:
            continue
        for u in G.predecessors(x):  # u->x in original
            Eb = G.edges[u, x].get('Eb', float('inf'))
            nd = d + Eb
            if nd < dist[u]:
                dist[u] = nd
                prev[u] = x
                heapq.heappush(pq, (nd, u))

    parent = {'D': None}
    for n in G.nodes():
        if n == 'D':
            continue
        parent[n] = prev[n]  # may be None if unreachable
    return parent

def ensure_feasible_parent(G, parent, sources):
    """Patch any None parents to the best feasible successor to avoid broken routes."""
    for s in sources:
        if parent.get(s) is None:
            feas = [(G.edges[s, v]['Eb'], v) for v in G.successors(s)
                    if G.edges[s, v]['Eb'] < float('inf') and not creates_cycle(parent, s, v)]
            if feas:
                feas.sort()
                parent[s] = feas[0][1]
    return parent

def build_k_hop_neighbors(G, k):
    """Return dict: node -> set of nodes within k hops (undirected reachability)."""
    UG = G.to_undirected()
    out = {}
    for s in G.nodes():
        vis = {s}
        q = deque([(s, 0)])
        while q:
            u, d = q.popleft()
            if d == k:
                continue
            for v in UG.neighbors(u):
                if v not in vis:
                    vis.add(v)
                    q.append((v, d+1))
        vis.discard(s)
        out[s] = vis
    return out

def best_response(G, pos, sources, parent, W_raw_full, mode="CAR",
                  alpha_f=0.5, corr_scale=3.0, B=None, neighbors=None, W_cache=None):
    """
    One sweep over players. Try feasible next-hops and apply strictly improving move.
    Uses caching to avoid recomputing W_agg excessively.
    Returns (improved_any, W_cache_out).
    """
    improved_any = False
    if W_cache is None:
        W_cache, _ = compute_aggregated_rates(G, pos, parent, W_raw_full, alpha_f, corr_scale)

    for i in sorted(sources):
        # candidate next-hops
        cand = set(G.successors(i))
        if neighbors is not None:
            cand &= set(neighbors.get(i, []))
        candidates = []
        for v in cand:
            if v != i and not creates_cycle(parent, i, v) and G.edges[i, v].get('Eb', float('inf')) < float('inf'):
                candidates.append(v)
        if not candidates:
            continue

        routes, _ = build_routes_and_users(parent)
        if mode == "MER":
            cur_cost = player_cost_MER(G, parent, routes, i)
        else:
            # Build pruned W_raw once for i
            descendants = []
            stack = [i]
            while stack:
                u = stack.pop()
                descendants.append(u)
                for v2, p2 in parent.items():
                    if p2 == u:
                        stack.append(v2)
            W_raw_wo = dict(W_raw_full)
            for u in descendants:
                if u != 'D':
                    W_raw_wo[u] = 0.0
            # current and pruned W_agg (once)
            W_with = W_cache
            W_without, _ = compute_aggregated_rates(G, pos, parent, W_raw_wo, alpha_f, corr_scale)
            route_nodes = path_to_root(parent, i)
            cur_cost = player_cost_CAR_from_cached(G, parent, route_nodes, W_with, W_without)

        best_v = None
        best_cost = cur_cost
        best_W = None

        for v in candidates:
            old_p = parent[i]
            parent[i] = v

            # Recompute W_agg *once* for tentative parent (shared for capacity + eval)
            W_agg_tmp, _ = compute_aggregated_rates(G, pos, parent, W_raw_full, alpha_f, corr_scale)
            if not capacity_ok(G, parent, W_agg_tmp, B=B):
                parent[i] = old_p
                continue

            if mode == "MER":
                routes_tmp, _ = build_routes_and_users(parent)
                new_cost = player_cost_MER(G, parent, routes_tmp, i)
            else:
                # Reuse previously computed W_raw_wo (descendants didn’t change)
                W_without_tmp, _ = compute_aggregated_rates(G, pos, parent, W_raw_wo, alpha_f, corr_scale)
                route_nodes = path_to_root(parent, i)
                new_cost = player_cost_CAR_from_cached(G, parent, route_nodes, W_agg_tmp, W_without_tmp)

            if new_cost + 1e-12 < best_cost:
                best_cost = new_cost
                best_v = v
                best_W = W_agg_tmp  # keep the computed W_agg for accepted move

            parent[i] = old_p

        if best_v is not None:
            parent[i] = best_v
            W_cache = best_W if best_W is not None else W_cache
            improved_any = True

    return improved_any, W_cache

def run_game_dynamics(G, pos, sources,
                      W_raw_val=200.0, Pc_target=0.99,
                      alpha_f=0.5, corr_scale=3.0,
                      B=None, mode="CAR", eps=1e-9, max_iters=200, k_hop=None):
    """
    Generic best-response driver. mode in {"MER","CAR"}.
    Returns (parent, W_agg, potential_value).
    """
    compute_link_metrics_tdma(G, pos, Pc_target=Pc_target)

    # W_raw: only sources generate raw rate
    W_raw = {n: 0.0 for n in G.nodes()}
    for s in sources:
        W_raw[s] = W_raw_val
    W_raw['D'] = 0.0

    # Initialization: MER (paper's recommendation) + feasibility patch
    parent = dijkstra_mer_parent(G, sink='D')
    parent = ensure_feasible_parent(G, parent, sources)

    # Optional: restrict candidates to k-hop neighborhood
    neigh = build_k_hop_neighbors(G, k_hop) if k_hop is not None else None

    # Iterate best responses with potential-based stopping
    W_agg, _ = compute_aggregated_rates(G, pos, parent, W_raw, alpha_f, corr_scale)
    prev_P = potential_energy(G, parent, W_agg)
    W_cache = W_agg

    for _it in range(max_iters):
        changed, W_cache = best_response(G, pos, sources, parent, W_raw,
                                         mode=mode, alpha_f=alpha_f, corr_scale=corr_scale, B=B,
                                         neighbors=neigh, W_cache=W_cache)
        P = potential_energy(G, parent, W_cache)

        # Stop if no change or tiny potential drop
        if (not changed) or (prev_P - P <= eps):
            break
        prev_P = P

    return parent, W_cache, prev_P

# Performance Comparisons and Plots

def compare_algorithms(G, pos, sources, W_raw_val=200.0, Pc_target=0.99,
                      alpha_f=0.5, corr_scale=3.0, B=None, max_enumeration_sources=6):
    """
    Compare MER, CAR, and Brute Force algorithms in terms of execution time and energy efficiency.
    """
    results = {}

    print("=" * 60)
    print("ALGORITHM COMPARISON")
    print("=" * 60)

    # Check if we can run brute force (only for small networks)
    can_run_brute_force = len(sources) <= max_enumeration_sources

    # 1. Brute Force (if applicable)
    if can_run_brute_force:
        print(f"\n1. BRUTE FORCE (N={len(sources)} sources)")
        print("-" * 40)
        start_time = time.time()
        brute_results = enumerate_all_trees(
            G, pos, sources,
            W_raw_val=W_raw_val,
            Pc_target=Pc_target,
            alpha_f=alpha_f,
            corr_scale=corr_scale,
            B=B,
            verbose=False
        )
        brute_time = time.time() - start_time

        if brute_results:
            best_brute = brute_results[0]
            results['brute_force'] = {
                'energy': best_brute['energy'],
                'time': brute_time,
                'parent': best_brute['parent'],
                'routes': best_brute['routes'],
                'total_solutions': len(brute_results)
            }
            print(f"   Execution time: {brute_time:.4f} seconds")
            print(f"   Total energy: {best_brute['energy']:.6e} J/s")
            print(f"   Total valid solutions found: {len(brute_results)}")
        else:
            print("   No valid solutions found")
    else:
        print(f"\n1. BRUTE FORCE - SKIPPED (N={len(sources)} > {max_enumeration_sources})")
        print("   Brute force enumeration skipped for large networks")

    # 2. MER Algorithm
    print(f"\n2. MER (Minimum Energy Routing)")
    print("-" * 40)
    start_time = time.time()
    parent_mer, W_agg_mer, _ = run_game_dynamics(
        G, pos, sources,
        W_raw_val=W_raw_val, Pc_target=Pc_target,
        alpha_f=alpha_f, corr_scale=corr_scale,
        B=B, mode="MER", eps=1e-9, max_iters=200, k_hop=None
    )
    mer_time = time.time() - start_time
    E_mer = total_energy(G, parent_mer, W_agg_mer)
    results['mer'] = {
        'energy': E_mer,
        'time': mer_time,
        'parent': parent_mer,
        'routes': extract_routes(parent_mer)
    }
    print(f"   Execution time: {mer_time:.4f} seconds")
    print(f"   Total energy: {E_mer:.6e} J/s")

    # 3. CAR Algorithm
    print(f"\n3. CAR (Correlation-Aware Routing)")
    print("-" * 40)
    start_time = time.time()
    parent_car, W_agg_car, _ = run_game_dynamics(
        G, pos, sources,
        W_raw_val=W_raw_val, Pc_target=Pc_target,
        alpha_f=alpha_f, corr_scale=corr_scale,
        B=B, mode="CAR", eps=1e-9, max_iters=200, k_hop=None
    )
    car_time = time.time() - start_time
    E_car = total_energy(G, parent_car, W_agg_car)
    results['car'] = {
        'energy': E_car,
        'time': car_time,
        'parent': parent_car,
        'routes': extract_routes(parent_car)
    }
    print(f"   Execution time: {car_time:.4f} seconds")
    print(f"   Total energy: {E_car:.6e} J/s")

    # Summary comparison
    print(f"\n4. COMPARISON SUMMARY")
    print("-" * 40)

    if can_run_brute_force and 'brute_force' in results:
        print(f"   Brute Force: {results['brute_force']['energy']:.6e} J/s ({results['brute_force']['time']:.4f}s)")
        print(f"   MER:         {results['mer']['energy']:.6e} J/s ({results['mer']['time']:.4f}s)")
        print(f"   CAR:         {results['car']['energy']:.6e} J/s ({results['car']['time']:.4f}s)")

        optimal_energy = results['brute_force']['energy']
        mer_gap = ((results['mer']['energy'] - optimal_energy) / optimal_energy) * 100
        car_gap = ((results['car']['energy'] - optimal_energy) / optimal_energy) * 100
        print(f"\n   Optimality gaps (vs Brute Force):")
        print(f"   MER: {mer_gap:.2f}%")
        print(f"   CAR: {car_gap:.2f}%")

        mer_speedup = results['brute_force']['time'] / results['mer']['time']
        car_speedup = results['brute_force']['time'] / results['car']['time']
        print(f"\n   Speedup factors (vs Brute Force):")
        print(f"   MER: {mer_speedup:.1f}x faster")
        print(f"   CAR: {car_speedup:.1f}x faster")
    else:
        print(f"   MER: {results['mer']['energy']:.6e} J/s ({results['mer']['time']:.4f}s)")
        print(f"   CAR: {results['car']['energy']:.6e} J/s ({results['car']['time']:.4f}s)")
        if results['mer']['energy'] < results['car']['energy']:
            car_gap = ((results['car']['energy'] - results['mer']['energy']) / results['mer']['energy']) * 100
            print(f"\n   CAR is {car_gap:.2f}% worse than MER")
        else:
            mer_gap = ((results['mer']['energy'] - results['car']['energy']) / results['car']['energy']) * 100
            print(f"\n   MER is {mer_gap:.2f}% worse than CAR")

    return results

def plot_comparison_results(results, G, pos):
    """
    Create comprehensive plots comparing the algorithms.
    """
    if not results:
        print("No results to plot")
        return

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Algorithm Comparison: MER vs CAR vs Brute Force', fontsize=16, fontweight='bold')

    # 1. Energy Comparison
    ax1 = axes[0, 0]
    algorithms = []
    energies = []
    colors = []

    if 'brute_force' in results:
        algorithms.append('Brute Force')
        energies.append(results['brute_force']['energy'])
        colors.append('green')

    algorithms.append('MER')
    energies.append(results['mer']['energy'])
    colors.append('blue')

    algorithms.append('CAR')
    energies.append(results['car']['energy'])
    colors.append('red')

    bars = ax1.bar(algorithms, energies, color=colors, alpha=0.7)
    ax1.set_ylabel('Total Energy (J/s)')
    ax1.set_title('Energy Efficiency Comparison')
    ax1.tick_params(axis='x', rotation=45)
    for bar, energy in zip(bars, energies):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height, f'{energy:.2e}', ha='center', va='bottom', fontsize=9)

    # 2. Execution Time Comparison
    ax2 = axes[0, 1]
    times = []
    time_algorithms = []
    time_colors = []

    if 'brute_force' in results:
        time_algorithms.append('Brute Force')
        times.append(results['brute_force']['time'])
        time_colors.append('green')

    time_algorithms.append('MER')
    times.append(results['mer']['time'])
    time_colors.append('blue')

    time_algorithms.append('CAR')
    times.append(results['car']['time'])
    time_colors.append('red')

    bars2 = ax2.bar(time_algorithms, times, color=time_colors, alpha=0.7)
    ax2.set_ylabel('Execution Time (seconds)')
    ax2.set_title('Execution Time Comparison')
    ax2.tick_params(axis='x', rotation=45)
    for bar, time_val in zip(bars2, times):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height, f'{time_val:.4f}s', ha='center', va='bottom', fontsize=9)

    # 3. MER Routing Tree
    ax3 = axes[1, 0]
    plot_tree_on_axis(G, pos, results['mer']['parent'], ax3,
                      f"MER Routing Tree\n(Energy: {results['mer']['energy']:.2e} J/s)")

    # 4. CAR Routing Tree
    ax4 = axes[1, 1]
    plot_tree_on_axis(G, pos, results['car']['parent'], ax4,
                      f"CAR Routing Tree\n(Energy: {results['car']['energy']:.2e} J/s)")

    plt.tight_layout()
    plt.show()

def run_scalability_analysis(max_sources=8, sink_neighbors=2):
    """
    Run scalability analysis to see how algorithms perform with different network sizes.
    """
    print("=" * 60)
    print("SCALABILITY ANALYSIS")
    print("=" * 60)

    scalability_data = {
        'sources': [],
        'mer_times': [],
        'car_times': [],
        'brute_times': [],
        'mer_energies': [],
        'car_energies': [],
        'brute_energies': []
    }

    for num_sources in range(3, max_sources + 1):
        print(f"\nTesting with {num_sources} sources...")

        # Build network
        G, pos, sources = pentagon(
            sink_pos=(5, 7), pentagon_center=(5, 5), radius=3,
            num_sources=num_sources, sink_neighbors=sink_neighbors
        )

        # Run comparison
        results = compare_algorithms(G, pos, sources, max_enumeration_sources=6)

        # Store data
        scalability_data['sources'].append(num_sources)
        scalability_data['mer_times'].append(results['mer']['time'])
        scalability_data['car_times'].append(results['car']['time'])
        scalability_data['mer_energies'].append(results['mer']['energy'])
        scalability_data['car_energies'].append(results['car']['energy'])

        if 'brute_force' in results:
            scalability_data['brute_times'].append(results['brute_force']['time'])
            scalability_data['brute_energies'].append(results['brute_force']['energy'])
        else:
            scalability_data['brute_times'].append(None)
            scalability_data['brute_energies'].append(None)

    plot_scalability_results(scalability_data)
    return scalability_data

def plot_scalability_results(data):
    """
    Plot scalability analysis results.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Time scalability
    ax1.plot(data['sources'], data['mer_times'], 'bo-', label='MER', linewidth=2, markersize=8)
    ax1.plot(data['sources'], data['car_times'], 'ro-', label='CAR', linewidth=2, markersize=8)

    # Filter out None values for brute force
    brute_sources = [s for s, t in zip(data['sources'], data['brute_times']) if t is not None]
    brute_times = [t for t in data['brute_times'] if t is not None]
    if brute_times:
        ax1.plot(brute_sources, brute_times, 'go-', label='Brute Force', linewidth=2, markersize=8)

    ax1.set_xlabel('Number of Sources')
    ax1.set_ylabel('Execution Time (seconds)')
    ax1.set_title('Execution Time Scalability')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Energy scalability
    ax2.plot(data['sources'], data['mer_energies'], 'bo-', label='MER', linewidth=2, markersize=8)
    ax2.plot(data['sources'], data['car_energies'], 'ro-', label='CAR', linewidth=2, markersize=8)

    brute_energies = [e for e in data['brute_energies'] if e is not None]
    if brute_energies:
        ax2.plot(brute_sources, brute_energies, 'go-', label='Brute Force', linewidth=2, markersize=8)

    ax2.set_xlabel('Number of Sources')
    ax2.set_ylabel('Total Energy (J/s)')
    ax2.set_title('Energy Efficiency Scalability')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Build network
    G, pos, sources = pentagon(
        sink_pos=(5, 7), pentagon_center=(5, 5), radius=3, num_sources=5, sink_neighbors=2
    )

    # Run comprehensive comparison
    print("Running comprehensive algorithm comparison...")
    results = compare_algorithms(G, pos, sources)

    # Plot comparison results
    plot_comparison_results(results, G, pos)

    # Run scalability analysis 
    print("\n" + "="*60)
    print("Running scalability analysis...")
    scalability_data = run_scalability_analysis(max_sources=7)