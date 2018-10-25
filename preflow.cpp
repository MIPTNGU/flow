#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <ctime>
#include <list>

typedef int TVector;
typedef int TFlow;

struct NetworkEdge
{
	//typedef int TVector;
	TVector start, finish;
	TFlow capacity, flow;
	bool full()
	{
		return capacity == flow;
	}
	bool can_push(TFlow flow_)
	{
		return flow + flow_ <= capacity;
	}
	void push(TFlow flow_)
	{
		if (can_push(flow_))
		{
			flow += flow_;
		}
	}
	bool equ(const NetworkEdge& other)
	{
		return start == other.start && finish == other.finish;
	}
	bool operator<(const NetworkEdge& other)
	{
		return std::make_pair(start, finish) < std::make_pair(other.start, other.finish);
	}
};

class Network
{
private:
	int vertices_number;
	TVector source, sink;
	std::vector < NetworkEdge > edges;
	std::vector < int > last_edge;
	std::vector < int > prev_edge;
	enum Empty_Edge
	{
		empty = -1
	};
public:
	~Network() = default;
	Network(size_t vertices_number_, size_t source_, size_t sink_,
		std::vector<NetworkEdge> edges_ = std::vector<NetworkEdge>()) :
		vertices_number(vertices_number_), source(source_), sink(sink_)
	{
		edges.clear();
		last_edge.resize(vertices_number);
		fill(last_edge.begin(), last_edge.end(), empty);
		prev_edge.clear(); 
		sort(edges_.begin(), edges_.end());
		for (NetworkEdge edge: edges_)
		{
			bool result = false;
			for (int edge_number = 0; edge_number < edges.size(); edge_number++)
			{
				if (edge.equ(edges[edge_number]))
				{
					edges[edge_number].capacity += edge.capacity;
					result = true;
					break;
				}
			}
			if (!result)
			{
				addEdge(edge.start, edge.finish, edge.capacity);
			}
		}
	}
	void updateEdge(int number, TFlow capacity)
	{
		edges[number].capacity += capacity;
	}
	void addEdge(int start, int finish, TFlow capacity)
	{
		edges.push_back({ start, finish, capacity, 0 });
		edges.push_back({ finish, start, 0, 0 });
		prev_edge.push_back(last_edge[start]);
		last_edge[start] = edges.size() - 2;
		prev_edge.push_back(last_edge[finish]);
		last_edge[finish] = edges.size() - 1;
	}
	int get_vertices_number()
	{
		return vertices_number;
	}
	int get_source()
	{
		return source;
	}
	int get_sink()
	{
		return sink;
	}
	int get_last(int vertex)
	{
		return last_edge[vertex];
	}
	std::vector < int > get_last_edges()
	{
		return last_edge;
	}
	TFlow get_flow()
	{
		TFlow current_flow = 0;
		for (int number_edge = last_edge[source]; number_edge != -1; number_edge = prev_edge[number_edge])
		{
			current_flow += edges[number_edge].flow;
		}
		return current_flow;
	}
	//typedef int iterator_type;
	friend class IOutputEdgeIterator;
	/*class IOutputEdgeIterator
	{	
		Network* network;
		iterator_type number_edge;
	public:
		IOutputEdgeIterator(Network* network_, iterator_type number_edge_) :
			network(network_), number_edge(number_edge_)
		{
		}
		bool valid()
		{
			return network != nullptr && number_edge != empty;
		}
		void next()
		{
			number_edge = network->prev_edge[number_edge];
		}
		TVector getStart(TVector vertex)
		{
			return network->edges[number_edge].start;
			//number_edge = last_edge[vertex];
		}
		TVector getFinish()
		{
			return network->edges[number_edge].finish;
		}
		TFlow getCapacity()
		{
			return network->edges[number_edge].capacity;
		}
		TFlow getFlow()
		{
			return network->edges[number_edge].flow;
		}
		//getResidual();
		NetworkEdge getEdge()
		{
			return network->edges[number_edge];
		}
		bool canPushFlow(TFlow flow)
		{
			return network->edges[number_edge].capacity - network->edges[number_edge].flow >= flow;
		}
		void pushFlow(TFlow flow)
		{
			network->edges[number_edge].flow += flow;
			network->edges[number_edge & 1].flow -= flow;
		}
		//for (it = begin(); it.valid(); it.next())
	};
	*/
};

typedef int iterator_type;
//friend class IOutputEdgeIterator;
class IOutputEdgeIterator
{
	Network* network;
	iterator_type number_edge;
public:
	IOutputEdgeIterator() : network(nullptr), number_edge(-1)
	{

	}
	IOutputEdgeIterator(Network* network_, iterator_type number_edge_) :
		network(network_), number_edge(number_edge_)
	{
	}
	//IOutputEdgeIterator& operator=(IOutputEdgeIterator& other)
	bool valid()
	{
		return network != nullptr && number_edge != -1;
	}
	void next()
	{
		number_edge = network->prev_edge[number_edge];
	}
	TVector get_start()
	{
		return network->edges[number_edge].start;
		//number_edge = last_edge[vertex];
	}
	int get_number_edge()
	{
		return number_edge;
	}
	TVector get_finish()
	{
		return network->edges[number_edge].finish;
	}
	TFlow get_capacity()
	{
		return network->edges[number_edge].capacity;
	}
	TFlow get_flow()
	{
		return network->edges[number_edge].flow;
	}
	//getResidual();
	NetworkEdge get_edge()
	{
		return network->edges[number_edge];
	}
	bool can_push_flow(TFlow flow)
	{
		return network->edges[number_edge].capacity - network->edges[number_edge].flow >= flow;
	}
	bool saturated()
	{
		return network->edges[number_edge].flow == network->edges[number_edge].capacity;
	}
	void push_flow(TFlow flow)
	{
		network->edges[number_edge].flow += flow;
		network->edges[number_edge ^ 1].flow -= flow;
	}
	//for (it = begin(); it.valid(); it.next())
};

class IAlgorithm
{
private:
	TFlow INF = 1e9;
	Network* graph;
	TFlow maxflow;
	TFlow path()
	{
		std::queue < int > q;
		q.push(graph->get_source());
		std::vector < int > num(graph->get_vertices_number());
		fill(num.begin(), num.end(), -1);
		std::vector < bool > used(graph->get_vertices_number());
		//par[s] = 0;
		//par1[s] = 0;
		//memset(used, false, graph.size() + 100);
		//memset(par, 0, graph.size() + 100);
		//memset(par1, 0, graph.size() + 100);
			used[graph->get_source()] = true;
			while (!q.empty())
			{
				int v = q.front();
				q.pop();
				if (v == graph->get_sink())
				{
					int u = graph->get_sink();
					TFlow f = INF;
					while (u != graph->get_source())
					{
						int number_edge = num[u];
						IOutputEdgeIterator it(graph, number_edge);
						f = std::min(f, it.get_capacity() - it.get_flow());
						//int num = par1[u];
						//f = std::min(f, graph[u1][num].c - graph[u1][num].f);
						u = it.get_start();
					}
					u = graph->get_sink();
					while (u != graph->get_source())
					{
						int number_edge = num[u];
						IOutputEdgeIterator it(graph, number_edge);
						it.push_flow(f);
						u = it.get_start();
					}
					return f;
				}
				for (IOutputEdgeIterator it(graph, graph->get_last(v)); it.valid(); it.next())
				{
					if (it.get_capacity() - it.get_flow() > 0 && !used[it.get_finish()])
					{
						used[it.get_finish()] = true;
						q.push(it.get_finish());
						num[it.get_finish()] = it.get_number_edge();
					}
				}
				/*for (int i = 0; i < graph[v].size(); i++)
				{
					if (graph[v][i].c - graph[v][i].f > 0 && !used[graph[v][i].u])
					{
						used[graph[v][i].u] = true;
						q.push(graph[v][i].u);
						par[graph[v][i].u] = v;
						par1[graph[v][i].u] = i;
					}
				}
				*/
			}
			return 0;
		}
public:
	IAlgorithm(Network* graph_) : graph(graph_), maxflow(0)
	{
	}
	~IAlgorithm() = default;
	TFlow run()
	{
		TFlow flow = path() + maxflow;
		if (flow == maxflow)
		{
			return maxflow;
		}
		else
		{
			maxflow = flow;
		}
		return run();
	}
};

class AlgorithmPreflow
{
private:
	TFlow INF = 1e9;
	Network* graph;
	TFlow maxflow;
	std::vector < TFlow > overflow;
	std::vector < int > height;
	std::vector < IOutputEdgeIterator > current_edge;
	//std::queue<int> q;
public:
	AlgorithmPreflow(Network* graph_) : graph(graph_), maxflow(0)
	{
		overflow.resize(graph->get_vertices_number());
		fill(overflow.begin(), overflow.end(), 0);
		height.resize(graph->get_vertices_number());
		fill(height.begin(), height.end(), 0);
		current_edge.resize(graph->get_vertices_number());
		for (int i = 0; i < graph->get_vertices_number(); i++)
		{
			current_edge[i] = IOutputEdgeIterator(graph, graph->get_last(i));
		}
		height[graph->get_source()] = graph->get_vertices_number();
	}
	~AlgorithmPreflow() = default;
	bool notfull(int v)
	{
		return overflow[v] > 0;
	}
	void push_flow(TFlow flow, IOutputEdgeIterator it)
	{
		it.push_flow(flow);
		overflow[it.get_start()] -= flow;
		overflow[it.get_finish()] += flow;
	}
	bool push(IOutputEdgeIterator& it)
	{
		int u = it.get_start();
		int v = it.get_finish();
		if (overflow[u] > 0 && height[u] == height[v] + 1 && !it.saturated())
		{
			TFlow flow = std::min(it.get_capacity() - it.get_flow(), overflow[u]);
			push_flow(flow, it);
			return true;
		}
		return false;
	}
	bool relabel(int u)
	{
		if (overflow[u] <= 0)
		{
			return false;
		}
		int minheight = INF;
		for (IOutputEdgeIterator it(graph, graph->get_last(u)); it.valid(); it.next())
		{
			int v = it.get_finish();
			if (!it.saturated())
			{
				if (height[v] < height[u])
				{
					return false;
				}
				else
				{
					minheight = std::min(minheight, height[v]);
				}
			}
		}
		height[u] = minheight + 1;
		return true;
	}
	void discharge(int u)
	{
		while (overflow[u] > 0)
		{
			if (!current_edge[u].valid())
			{
				relabel(u);
				current_edge[u] = { graph, graph->get_last(u) };
			}
			else
			{
				if (!current_edge[u].saturated() && height[u] == height[current_edge[u].get_finish()] + 1)
				{
					push(current_edge[u]);
				}
				else
				{
					current_edge[u].next();
				}
			}
		}
	}
	TFlow run()
	{
		int cur_flow = 0;
		for (IOutputEdgeIterator it(graph, graph->get_last(graph->get_source())); it.valid(); it.next())
		{
			overflow[graph->get_source()] += it.get_capacity();
			push_flow(it.get_capacity(), it);
		}
		while (true)
		{
			bool can_do_push_or_relabel = false;
			for (int i = 0; i < graph->get_vertices_number(); i++)
			{
				if (i != graph->get_source() && i != graph->get_sink() && overflow[i] > 0)
				{
					discharge(i);
					can_do_push_or_relabel = true;
				}
			}
			if (!can_do_push_or_relabel)
			{
				break;
			}
		}
		for (IOutputEdgeIterator it(graph, graph->get_last(graph->get_source())); it.valid(); it.next())
		{
			cur_flow += it.get_flow();
		}
		maxflow = cur_flow;
		return cur_flow;
	}
};

/*void bfs(int v, Network* graph)
{
	std::vector < bool > used(graph->get_vertices_number());
	std::queue < int > q;
	q.push(v);
	while (!q.empty())
	{
		int u = q.front;
		for (IOutputEdgeIterator it(graph, graph->get_last()); it.valid(); it.next())
		{
			if (!used)
		}
	}
}
*/

/*class LayeredNetwork
{
	Network graph;
	LayeredNetwork(Network* other) :
	{
		int source = other->get_source();
		int sink = other->get_sink();
		bfs(source);
	}
};
*/

class AlgorithmMalhotra
{
	Network* graph;
	TFlow maxflow;
	AlgorithmMalhotra(Network* graph_) : graph(graph_), maxflow(0)
	{

	}
};

int main()
{
#ifdef _DEBUG
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif
	while (true)
	{
		int n, s, t, m;
		std::cin >> n;
		if (n == 0)
		{
			return 0;
		}
		std::cin >> s >> t >> m;
		s--;
		t--;
		std::vector < NetworkEdge > pok;
		for (int i = 1; i <= m; i++)
		{
			int v, u, c;
			std::cin >> v >> u >> c;
			v--;
			u--;
			pok.push_back({ v, u, c, 0 });
			pok.push_back({ u, v, c, 0 });
		}
		Network graph(n, s, t, pok);
		AlgorithmPreflow mypok(&graph);
		std::cout << mypok.run() << "\n";
	}
	return 0;
}
