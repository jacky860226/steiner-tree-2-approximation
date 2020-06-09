# steiner-tree-2-approximation
steiner tree 2-approximation O(E + V log V)

Paper: [A faster approximation algorithm for the Steiner problem in graphs](https://www.sciencedirect.com/science/article/pii/002001908890066X)
## Example usage
```cpp
#include "Structure.hpp"
#include "Solver.hpp"
#include <iostream>
#include <random>
#include <vector>
using namespace std;
mt19937 MT(7122);
steiner_tree::UndirectedGraph<double> G;
int main()
{
    int N = 100, M = 10;
    vector<size_t> terminals;
    for (int i = 0; i < N; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            size_t eid = G.addEdge(i, j, MT()*0.001);
        }
        terminals.emplace_back(i);
    }
    shuffle(terminals.begin(), terminals.end(), MT);
    terminals.erase(terminals.begin() + M, terminals.end());
    steiner_tree::Solver<double> solver(G);
    auto res = solver.solve(terminals);
    for (auto eid : res)
    {
        cout << G.getEdge(eid).v1 << ' ' << G.getEdge(eid).v2 << ' ' << G.getEdge(eid).cost << endl;
    }
    return 0;
}
```