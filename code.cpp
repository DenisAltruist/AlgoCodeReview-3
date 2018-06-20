#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <utility>
#include <queue>
#include <cstdio>


struct Edge {
    int v, u, cap;
    long long cost;
    int flow, idx;
};

struct FlowDecomposition {
    int cost, flow;
    std::vector<std::vector<Edge>> pathWays;

    void addPath(const std::vector<Edge> path) {
        pathWays.push_back(path);
    }

    FlowDecomposition(): cost(0), flow(0){};
};

class Graph {
 private:
    const int INF = std::numeric_limits<int>::max();
    std::vector<Edge> edges;
    std::vector<std::vector<int>> g;
    int n;

    bool isSaturated(int idx) const {
        if (edges[idx].flow == edges[idx].cap) {
            return true;
        }
        return false;
    }

    bool isReversed(int idx) const {
        return (idx % 2 == 1);
    }

    Edge getEdge(int idx) {
        return edges[idx];
    }

    std::vector<int> getShortestPath(int s, int t) {
        auto res = dijkstra(s);
        std::vector<int> prev = res.second;
        if (prev[t] == -1) {
            return std::vector<int>();
        }
        std::vector<int> path;
        int curv = t;
        while (curv != s) {
            path.push_back(prev[curv]);
            curv = getNeighbor(curv, prev[curv]);
        }
        std::reverse(path.begin(), path.end());
        return path;
    }

    std::vector<int> getDistances(int s) {
        auto res = dijkstra(s);
        return res.first;
    }

    std::pair<std::vector<int>, std::vector<int>> dijkstra(int s) const {
        std::priority_queue<std::pair<int, int>> q;
        std::vector<int> dist(n + 1, INF);
        std::vector<int> prev(n + 1, -1);
        dist[s] = 0;
        for (int i = 1; i <= n; ++i) {
            q.push(std::make_pair(-dist[i], i));
        }
        while (!q.empty()) {
            int v = q.top().second;
            int curDist = -q.top().first;
            q.pop();
            if (curDist > dist[v] || dist[v] == INF) {
                continue;
            }
            for (int idx: g[v]) {
                if (isSaturated(idx)) {
                    continue;
                }
                int u = getNeighbor(v, idx);
                int costOfEdge = edges[idx].cost;
                if (dist[u] > dist[v] + costOfEdge) {
                    dist[u] = dist[v] + costOfEdge;
                    prev[u] = idx;
                    q.push(std::make_pair(-dist[u], u));
                }
            }
        }
        return std::make_pair(dist, prev);
    }

    std::vector<int> calcPotentials(int s) {
        std::vector<int> p(n + 1, INF);
        p[s] = 0;
        bool isSomethingRelaxed;
        do {
            isSomethingRelaxed = false;
            for (int v = 1; v <= n; ++v) {
                if (p[v] != INF) {
                    for (int idx: g[v]) {
                        if (isSaturated(idx)) {
                            continue;
                        }
                        int u = getNeighbor(v, idx);
                        if (p[u] > p[v] + edges[idx].cost) {
                            p[u] = p[v] + edges[idx].cost;
                            isSomethingRelaxed = true;
                        }
                    }
                }
            }

        } while (isSomethingRelaxed);
        return p;
    }

    void addEdge(const Edge& e) {
        g[e.v].push_back(edges.size());
        edges.push_back(e);
    }

    int getNeighbor(int v, int idx) const {
        if (edges[idx].v == v) {
            return edges[idx].u;
        }
        return edges[idx].v;
    }

    void pushFlow(int idx, int value) {
        if (idx % 2 == 0) {
            edges[idx].flow += value;
            edges[idx + 1].flow -= value;
        } else {
            edges[idx].flow += value;
            edges[idx - 1].flow -= value;
        }
    }

 public:
    Graph(int nValue, const std::vector<Edge>& edgesValue) {
        n = nValue;
        g.resize(n + 1);
        for (Edge curEdge: edgesValue) {
            addEdge(curEdge);
            addEdge(Edge{curEdge.u, curEdge.v, 0, -curEdge.cost, curEdge.flow, curEdge.idx});
        }
    }

    void reweight(const std::vector<int>& p) {
        for (Edge e: edges) {
            e.cost = e.cost + p[e.v] - p[e.u];
        }
    }

    int calcCost() const {
        int res = 0;
        for (Edge e: edges) {
            if (e.flow > 0) {
                res += e.cost;
            }
        }
        return res;
    }

    void getPathFlow(int v, std::vector<int>& prev) {
        for (int idx: g[v]) {
            if (isReversed(idx)) {
                continue;
            }
            int u = getNeighbor(v, idx);
            if (prev[u] == -1 && isSaturated(idx)) {
                prev[u] = idx;
                getPathFlow(u, prev);
            }
        }
    }

    FlowDecomposition getMinCostMaxFlow(int cntOfIters) {
        std::vector<int> p = calcPotentials(1);
        std::vector<Edge> reweightedEdges;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (i % 2 == 0) {
                reweightedEdges.push_back(edges[i]);
            }
        }
        FlowDecomposition res;
        Graph G(n, reweightedEdges);
        G.reweight(p);
        for (int it = 1; it <= cntOfIters; ++it) {
            std::vector<int> shortestPath = G.getShortestPath(1, n);
            std::vector<int> dist = G.getDistances(1);
            if (shortestPath.empty()) {
                break;
            }
            for (int idx: shortestPath) {
                G.pushFlow(idx, 1);
            }
            res.flow = it;
            for (size_t i = 0; i < p.size(); ++i) {
                p[i] = -p[i];
            }
            G.reweight(p);
            for (int i = 1; i <= n; ++i) {
                if (dist[i] == INF) {
                    continue;
                }
                p[i] = -p[i] + dist[i];
            }
            G.reweight(p);
        }
        res.cost = G.calcCost();
        for (int j = 1; j <= res.flow; ++j) {
            std::vector<int> prev(n + 1, -1);
            G.getPathFlow(1, prev);
            int curv = n;
            std::vector<Edge> pathWay;
            while (curv != 1) {
                int idx = prev[curv];
                G.pushFlow(idx, -1);
                pathWay.push_back(G.getEdge(idx));
                curv = G.getNeighbor(curv, idx);
            }
            std::reverse(pathWay.begin(), pathWay.end());
            res.addPath(pathWay);
        }
        return res;
    }

};

int main() {
    std::ios_base::sync_with_stdio(false);
    //freopen("input.txt", "r", stdin);
    int n, m, k;
    std::cin >> n >> m >> k;
    std::vector<Edge> edges;
    for (int i = 1; i <= m; ++i) {
        int v, u, cost;
        std::cin >> v >> u >> cost;
        if (v == u) {
            continue;
        }
        edges.push_back(Edge{v, u, 1, cost, 0, i});
        edges.push_back(Edge{u, v, 1, cost, 0, i});
    }
    Graph G(n, edges);
    FlowDecomposition res = G.getMinCostMaxFlow(k);
    if (res.flow < k) {
        std::cout << "-1\n";
    } else {
        std::cout << std::fixed << std::setprecision(10) << (res.cost / (k + 0.0)) << "\n";
        for (size_t i = 0; i < res.pathWays.size(); ++i) {
            std::vector<Edge> curPath = res.pathWays[i];
            std::cout << curPath.size();
            for (auto e: curPath) {
                std::cout << " " << e.idx;
            }
            std::cout << "\n";
        }
    }
    return 0;
}
