//Ricardo Andres Arriaga Quezada
//A01570553

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <float.h>
#include <algorithm>
#include <climits>
#include <map>

using namespace std;

struct Colonia
{
    string name;
    float x, y;
    int c;
    
    Colonia()
    {
        name = "-";
        x = 0;
        y = 0;
        c = 0;
    }

    Colonia(string name, float x, float y, int c)
    {
        this->name = name;
        this->x = x;
        this->y = y;
        this->c = c;
    }    
};

struct Graph 
{
    int V, E, costKruskal = 0, costTSP = 0;
    map<string, int> colonias;
    vector<int> centrales;
    vector<int> noCentrales;
    vector<pair<int, pair<string, string>>> edges;
    vector<vector<int>> adjList;
    vector<vector<int>> colIntermedias;
    vector<pair<string, string>> kruskalAns;

    vector<Colonia> listaColonias;
    vector<Colonia> nuevasColonias;

    Graph(int V, int E)
    { 
        this->V = V; 
        this->E = E;
        adjList.resize(V, vector<int>(E, INT_MAX));
        colIntermedias.resize(V, vector<int>(E, -1));
    }

    void addEdge(string u, string v, int w) 
    { 
        edges.push_back({w, {u, v}});
        adjList[colonias[u]][colonias[v]] = w;
        adjList[colonias[v]][colonias[u]] = w;
    } 

    void loadAdjacencyList();
    void kruskal(ofstream &output);
    void floyd(ofstream &output);
    void TSP(ofstream &output);
    void nuevaColoniaCercana(ofstream &output);
};

struct DisjointSets 
{ 
    int *parent, *rnk; 
    int n; 
  
    DisjointSets(int n)
    {
        this->n = n; 
        parent = new int[n+1]; 
        rnk = new int[n+1]; 
        
        for (int i = 0; i <= n; i++)
        { 
            rnk[i] = 0; 
            parent[i] = i; 
        } 
    } 
  
    // Find the parent of a node 'u' 
    // Path Compression 
    int find(int u) 
    { 
        /* Make the parent of the nodes in the path 
           from u--> parent[u] point to parent[u] */
        if (u != parent[u])
        {
            parent[u] = find(parent[u]); 
        }
        return parent[u]; 
    } 
  
    // Union by rank 
    void merge(int x, int y) 
    { 
        x = find(x), y = find(y); 
  
        /* Make tree with smaller height 
           a subtree of the other tree  */
        if (rnk[x] > rnk[y])
        {
            parent[y] = x; 
        }
        else // If rnk[x] <= rnk[y] 
        {
            parent[x] = y; 
        } 
            
        if(rnk[x] == rnk[y])
        {
            rnk[y]++; 
        }
    } 
}; 

void TSPPath(vector<vector<int>> inter, int u, int v, vector<int> &path) 
{
    if(inter[u][v] != -1)
    {
        TSPPath(inter, u, inter[u][v], path);
        path.push_back(inter[u][v]);
    }
    else
    {
        path.push_back(u);
    }
}

void recrearPath(vector<vector<int>> inter, int u, int v, vector<int> &path)
{
    if(inter[u][v] != -1)
    {
        recrearPath(inter, u, inter[u][v], path);
        path.push_back(inter[u][v]);
        recrearPath(inter, inter[u][v], v, path);
    }
}

float distancia(Colonia &c1, Colonia &c2)
{
    return sqrt((c1.x - c2.x) * (c1.x - c2.x) + (c1.y - c2.y) * (c1.y - c2.y));
}

//Complejidad O(n^3)
void Graph::loadAdjacencyList()
{
    for(int i = 0; i < V; i++)
    {
        for(int j = 0 ; j < V; j++)
        {
            for(int k = 0; k < V; k++)
            {
                if(adjList[j][i] != INT_MAX && adjList[i][k] != INT_MAX && adjList[j][i] + adjList[i][k] < adjList[j][k])
                {
                    adjList[j][k] = adjList[j][i] + adjList[i][k];
                    colIntermedias[j][k] = i;
                }
            }
        }
    }
}

//Complejidad O(E Log2 E)
void Graph::kruskal(ofstream &output)
{ 
    output << "------------------------------" << endl;
    output << "1 - Cableado optimo de nueva conexion." << endl;

    sort(edges.begin(), edges.end());
    DisjointSets ds(V);

    for(auto it:edges)
    {
        int u = colonias[it.second.first];
        int v = colonias[it.second.second];
        int set_u = ds.find(u);
        int set_v = ds.find(v);

        if(set_u != set_v)
        {
            ds.merge(u, v);
            costKruskal += it.first;
            output << it.second.first << " - " << it.second.second << " : " << it.first << endl; 
        }
    }
    output << endl;
    output << "Costo Total: " << costKruskal << endl;
    output << "------------------------------" << endl;
}

void Graph::TSP(ofstream &output)
{
    output << "2 - La ruta optima." << endl;
    
    vector<int> pathTSP;
    string ans;

    for(int i = 0; i <= noCentrales.size() - 1; i++)
    {
        
        if(i+1 == noCentrales.size())
        {
            costTSP += adjList[noCentrales[i]][noCentrales[0]];
            TSPPath(colIntermedias, noCentrales[i], noCentrales[0], pathTSP);
        }
        else
        {
            costTSP += adjList[noCentrales[i]][noCentrales[i+1]];
            TSPPath(colIntermedias, noCentrales[i], noCentrales[i+1], pathTSP);
        }
    }
    
    for(int i = 0; i < pathTSP.size(); i++)
    {
        ans += listaColonias[pathTSP[i]].name + " - ";
    }

    ans += listaColonias[noCentrales[0]].name;

    output << ans << endl;
    output << "La Ruta Optima tiene un costo total de: " << costTSP << endl;
    output << "------------------------------" << endl;
}

void Graph::floyd(ofstream &output)
{
    string ans;
    vector<int> pathFloyd;

    output << "3 - Caminos mas cortos entre centrales." << endl;

    for(int i = 0; i < centrales.size() - 1; i++)
    {
        for(int j = i+1; j < centrales.size(); j++)
        {
            if(adjList[i][j] == INT_MAX)
            {
                output << "No Path" << endl;
            }
            else
            {
                pathFloyd.push_back(centrales[i]);
                recrearPath(colIntermedias, centrales[i], centrales[j], pathFloyd);
                
                for(int k = 0; k < pathFloyd.size(); k++)
                {
                    ans += listaColonias[pathFloyd[k]].name + " - ";
                }

                ans += listaColonias[centrales[j]].name;

                output << ans << " (" << adjList[centrales[i]][centrales[j]] << ") " << endl;
            }
            ans = "";
            pathFloyd.clear();
        }
    }
    output << "------------------------------" << endl;
}

//Complejidad O(n*m)
void Graph::nuevaColoniaCercana(ofstream &output)
{
    output << "4 - Conexion de nuevas colonias." << endl;

    for(int i = 0; i < nuevasColonias.size(); i++)
    {
        Colonia ans1;
        Colonia ans2;
        float minDist = FLT_MAX;

        for(int j = 0; j < listaColonias.size(); j++)
        {
            float dist = distancia(nuevasColonias[i], listaColonias[j]);

            if(dist < minDist)
            {
                minDist = dist;
                ans1 = nuevasColonias[i];
                ans2 = listaColonias[j];
            }
        }
        output << ans1.name << " debe conectarse con " << ans2.name << endl;
    }
    output << "------------------------------" << endl;
}

int main()
{
    int n, m, q, c, w;
    float x, y;
    string nombre, a, b;
    
    ifstream file("in02.txt");
    ofstream output;
    
    file >> n >> m >> q;

    Graph g(n, m);

    for(int i = 0; i < n; i++)
    {
        file >> nombre >> x >> y >> c;
        g.listaColonias.push_back(Colonia(nombre, x, y, c));
        g.colonias.insert(pair<string, int>(nombre, i));
        
        if(c == 1)
        {
            g.centrales.push_back(i);
        }
        else
        {
            g.noCentrales.push_back(i);
        }
    }

    for (int i = 0; i < n; i++)
    {
        g.adjList[i][i] = 0;
        for (int j = i+1; j < m; j++)
        {
            g.adjList[i][j] = INT_MAX;
        }
    }

    for (int i = 0; i < m; i++)
    {
        file >> a >> b >> w;
        g.addEdge(a, b, w);
    }

    for(int i = 0; i < q; i++)
    {
        file >> nombre >> x >> y;
        g.nuevasColonias.push_back(Colonia(nombre, x, y, 0));
    }

    file.close();

    output.open("checking2.txt");

    g.loadAdjacencyList();
    g.kruskal(output);
    g.TSP(output);
    g.floyd(output);
    g.nuevaColoniaCercana(output);

    output.close();

    return 0;
}