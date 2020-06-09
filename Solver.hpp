#pragma once
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <tuple>
#include <unordered_set>
#include <stack>
#include <cassert>

#include "Structure.hpp"
#include "Queue.hpp"
#include "DisjoinSet.hpp"
namespace steiner_tree
{
    template <class CostTy>
    class Solver
    {
        const size_t INVLID;
        const CostTy INF;

        const UndirectedGraph<CostTy> &G;

        std::vector<CostTy> dist;
        std::vector<size_t> prev_eid;
        std::vector<size_t> index;

        std::vector<size_t> &deg = prev_eid; // memory reduce

        std::vector<std::tuple<CostTy, size_t>> CrossEdge;
        std::vector<size_t> SelectKEdge;
        std::unordered_set<size_t> used_eid;

        DisjoinSet ds;

        void SPFA(const std::vector<size_t> &terminalVertice)
        {
            auto N = G.getVertexNum();
            std::vector<bool> inqueue(N, false);
            Queue<size_t> q;
            dist.reserve(N);
            prev_eid.reserve(N);
            index.reserve(N);
            for (size_t i = 0; i < N; ++i)
            {
                dist.emplace_back(INF);
                prev_eid.emplace_back(INVLID);
                index.emplace_back(INVLID);
            }
            for (auto v : terminalVertice)
            {
                dist[v] = 0;
                index[v] = v;
                q.push_back(v);
                inqueue[v] = true;
            }
            while (!q.empty())
            {
                auto v = q.front();
                q.pop_front();
                inqueue[v] = false;
                for (auto eid : G.getAdjacencyList(v))
                {
                    const auto &edge = G.getEdge(eid);
                    const auto dual = edge.getDual(v);
                    if (dist[dual] > dist[v] + edge.cost)
                    {
                        dist[dual] = dist[v] + edge.cost;
                        prev_eid[dual] = eid;
                        index[dual] = index[v];
                        if (!inqueue[dual])
                        {
                            if (!q.empty() && dist[dual] <= dist[q.front()])
                                q.push_front(dual);
                            else
                                q.push_back(dual);
                            inqueue[dual] = true;
                        }
                    }
                }
            }
        }

        void calculateCrossEdge()
        {
            size_t eid = 0;
            for (const auto &edge : G.getEdges())
            {
                if (index[edge.v1] != index[edge.v2])
                {
                    CrossEdge.emplace_back(dist[edge.v1] + dist[edge.v2] + edge.cost, eid);
                }
                eid++;
            }
            std::sort(CrossEdge.begin(), CrossEdge.end());
        }

        void Kruskal_1()
        {
            SelectKEdge.clear();
            ds.init(G.getVertexNum());
            for (const auto &TUS : CrossEdge)
            {
                size_t eid;
                std::tie(std::ignore, eid) = TUS;
                const auto &edge = G.getEdge(eid);
                if (!ds.same(index[edge.v1], index[edge.v2]))
                {
                    SelectKEdge.emplace_back(eid);
                    ds.Union(index[edge.v1], index[edge.v2]);
                }
            }
        }

        void edgeRecover()
        {
            std::vector<bool> used(G.getEdgeNum(), false);
            CrossEdge.clear();
            for (auto eid : SelectKEdge)
            {
                const auto &edge = G.getEdge(eid);
                CrossEdge.emplace_back(edge.cost, eid);
                for (int t = 0; t < 2; ++t)
                {
                    size_t v = (t & 1) ? edge.v1 : edge.v2;
                    while (v != index[v] && !used[prev_eid[v]])
                    {
                        used[prev_eid[v]] = true;
                        const auto &prev_edge = G.getEdge(prev_eid[v]);
                        CrossEdge.emplace_back(prev_edge.cost, prev_eid[v]);
                        v = prev_edge.getDual(v);
                    }
                }
            }
        }

        void Kruskal_2()
        {
            std::fill(deg.begin(), deg.end(), 0);
            used_eid.clear();
            std::sort(CrossEdge.begin(), CrossEdge.end());
            ds.init(G.getVertexNum());
            for (const auto &TUS : CrossEdge)
            {
                size_t eid;
                std::tie(std::ignore, eid) = TUS;
                const auto &edge = G.getEdge(eid);
                if (!ds.same(edge.v1, edge.v2))
                {
                    used_eid.insert(eid);
                    ++deg[edge.v1];
                    ++deg[edge.v2];
                    ds.Union(edge.v1, edge.v2);
                }
            }
        }

        bool edgeReduce(const std::vector<size_t> &terminalVertice)
        {
            size_t N = G.getVertexNum();
            std::stack<size_t> stack;
            std::vector<bool> is_terminal(N, false);
            for (auto v : terminalVertice)
            {
                is_terminal[v] = true;
            }
            bool isConnected = true;
            for (size_t i = 0; i < N; ++i)
            {
                if (is_terminal[i] && deg[i] == 0)
                    isConnected = false;
                if (!is_terminal[i] && deg[i] == 1)
                {
                    deg[i]--;
                    stack.emplace(i);
                }
            }
            if (terminalVertice.size() <= 1)
                isConnected = true;

            if (!isConnected)
                return false;

            while (!stack.empty())
            {
                auto v = stack.top();
                stack.pop();

                for (auto eid : G.getAdjacencyList(v))
                {
                    if (used_eid.find(eid) == used_eid.end())
                        continue;
                    used_eid.erase(eid);
                    auto next = G.getEdge(eid).getDual(v);
                    deg[next]--;
                    if (deg[next] == 0 && !is_terminal[next])
                        stack.emplace(next);
                }
            }
            SelectKEdge.clear();
            std::copy(used_eid.begin(), used_eid.end(), std::back_inserter(SelectKEdge));
        }

    public:
        Solver(const UndirectedGraph<CostTy> &G, const CostTy INF = std::numeric_limits<CostTy>::max() / 2 - 1) : G(G), INVLID(std::numeric_limits<size_t>::max()), INF(INF)
        {
            assert(G.getEdgeCosts() < INF);
        }
        std::vector<size_t> solve(const std::vector<size_t> &terminalVertice)
        {
            SPFA(terminalVertice);
            calculateCrossEdge();
            Kruskal_1();
            edgeRecover();
            Kruskal_2();
            edgeReduce(terminalVertice);
            return SelectKEdge;
        }
    };
} // namespace steiner_tree