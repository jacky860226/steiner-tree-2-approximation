#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <memory>
#include <stack>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "DisjoinSet.hpp"
#include "Queue.hpp"
#include "Structure.hpp"
namespace steiner_tree
{
template <class CostTy> class Solver
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

        dist.clear();
        dist.resize(N, INF);

        prev_eid.clear();
        prev_eid.resize(N, INVLID);

        index.clear();
        index.resize(N, INVLID);

        for (auto v : terminalVertice)
        {
            dist.at(v) = 0;
            index.at(v) = v;
            q.push_back(v);
            inqueue.at(v) = true;
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
        CrossEdge.clear();
        size_t eid = 0;
        for (const auto &edge : G.getEdges())
        {
            if (index[edge.v1] != index[edge.v2])
            {
                if (index[edge.v1] != INVLID && index[edge.v2] != INVLID)
                {
                    CrossEdge.emplace_back(dist[edge.v1] + dist[edge.v2] + edge.cost, eid);
                }
            }
            ++eid;
        }
    }

    bool Kruskal_1(size_t terminalSize)
    {
        SelectKEdge.clear();
        std::sort(CrossEdge.begin(), CrossEdge.end());
        ds.init(G.getVertexNum());
        size_t unionCnt = 0;
        for (const auto &TUS : CrossEdge)
        {
            size_t eid;
            std::tie(std::ignore, eid) = TUS;
            const auto &edge = G.getEdge(eid);
            if (!ds.same(index[edge.v1], index[edge.v2]))
            {
                SelectKEdge.emplace_back(eid);
                ds.Union(index[edge.v1], index[edge.v2]);
                ++unionCnt;
            }
        }
        return terminalSize == unionCnt + 1;
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

    void edgeReduce()
    {
        size_t N = G.getVertexNum();
        std::stack<size_t> stack;
        for (size_t i = 0; i < N; ++i)
        {
            if (index[i] != i && deg[i] == 1)
            {
                --deg[i];
                stack.emplace(i);
            }
        }

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
                --deg[next];
                if (deg[next] == 0 && index[next] != next)
                    stack.emplace(next);
            }
        }
    }

  public:
    Solver(const UndirectedGraph<CostTy> &G, const CostTy INF = std::numeric_limits<CostTy>::max() / 2 - 1)
        : INVLID(std::numeric_limits<size_t>::max()), INF(INF), G(G)
    {
        assert(G.getEdgeCosts() < INF);
    }
    std::shared_ptr<std::vector<size_t>> solve(const std::vector<size_t> &terminalVertice)
    {
        SPFA(terminalVertice);
        calculateCrossEdge();
        if (!Kruskal_1(terminalVertice.size()))
            return nullptr;
        edgeRecover();
        Kruskal_2();
        edgeReduce();
        SelectKEdge.clear();
        std::copy(used_eid.begin(), used_eid.end(), std::back_inserter(SelectKEdge));
        return std::make_shared<std::vector<size_t>>(SelectKEdge);
    }
};
} // namespace steiner_tree