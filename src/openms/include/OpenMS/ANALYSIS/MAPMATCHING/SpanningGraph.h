#pragma once
#include <stdio.h>
#include <stdlib.h>


using namespace OpenMS;
using namespace std;

class SpanningGraph
{  
    
public:
    int V;   
    vector<vector<int>> adj; 
    SpanningGraph(int V);   
    void addEdge(int v, int w);   
    bool containsCycle();  
    bool DFS(int v, int source, vector<bool>& visited); 
};
 
SpanningGraph::SpanningGraph(int V)
{
    this->V = V;
    vector<vector<int>> adj;
    for (Size i = 0; i < V; i++)
    {
      vector<int> a;
      adj.push_back(a);
    }
    this->adj = adj;
}
 
void SpanningGraph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}

bool SpanningGraph::containsCycle()
{
  vector<bool> visited(V);
  for (int i = 0; i < V; i++)
  {
    visited[i] = false;
  }
    
  for (int v = 0; v < V; v++)
  {
    if (!visited[v])
    {
      if (DFS(v, visited, -1))
      {
        return true;
      }
    }
  } 
  return false;

}

//DFS helper
bool SpanningGraph::DFS(int v, int source, vector<bool>& visited)
{
    visited[v] = true;
    for (unsigned int i = 0; i < adj[v].size(); i++)
    {
      if (!visited[adj[v][i]])
      {
        if (DFS(adj[v][i], visited, v)) {return true;}
      }
      else if (adj[v][i] != source) {return true;}
    }
    return false;
}

