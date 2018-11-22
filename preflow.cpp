#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <ctime>
#include <list>

typedef int TVertex;
typedef int TFlow;

struct NetworkEdge
{
	TVertex start, finish;
	TFlow capacity, flow;
	NetworkEdge(TVertex start_, TVertex finish_, TFlow capacity_, TFlow flow_ = 0) : start(start_), finish(finish_),
		capacity(capacity_), flow(flow_)
	{
	}
	~NetworkEdge() = default;
};

bool equal_ends(const NetworkEdge& first, const NetworkEdge& second)
{
	return first.start == second.start && first.finish == second.finish;
}

bool compare(const NetworkEdge& first, const NetworkEdge& second)
{
	return std::make_pair(first.start, first.finish) < std::make_pair(second.start, second.finish);
}



class Network
{
private:
	int vertices_number;
	TVertex source, sink;
	std::vector < NetworkEdge > edges;
	std::vector < int > last_edge;
	std::vector < int > prev_edge;
	std::vector < int > last_rev_edge;
	std::vector < int > prev_rev_edge;
	static const int empty = -1;
	void init()
	{
		edges.clear();
		last_edge.resize(vertices_number);
		fill(last_edge.begin(), last_edge.end(), empty);
		last_rev_edge.resize(vertices_number);
		fill(last_rev_edge.begin(), last_rev_edge.end(), empty);
		prev_edge.clear();
	}
public:
	~Network() = default;
	Network(size_t vertices_number_, size_t source_, size_t sink_,
		std::vector<NetworkEdge> edges_ = std::vector<NetworkEdge>()) :
		vertices_number(vertices_number_), source(source_), sink(sink_)
	{
		init();
		sort(edges_.begin(), edges_.end(), compare);
		for (NetworkEdge edge : edges_)
		{
			if (!edges.empty() && equal_ends(edge, edges.back()))
			{
				edges.back().capacity += edge.capacity;
			}
			else
			{
				addEdge(edge.start, edge.finish, edge.capacity);
			}
		}
	}
private:
	void updateEdge(int number, TFlow capacity)
	{
		edges[number].capacity += capacity;
	}
	void addEdge(TVertex start, TVertex finish, TFlow capacity)
	{
		add_edge_without_reverse(start, finish, capacity);
		add_edge_without_reverse(finish, start, 0);
	}
	void add_edge_without_reverse(TVertex start, TVertex finish, TFlow capacity)
	{
		edges.push_back({ start, finish, capacity });
		prev_edge.push_back(last_edge[start]);
		last_edge[start] = edges.size() - 1;
		prev_rev_edge.push_back(last_rev_edge[finish]);
		last_rev_edge[finish] = edges.size() - 1;
	}
public:
	int get_vertices_number() const
	{
		return vertices_number;
	}
	int get_edges_number() const
	{
		return edges.size();
	}
	TVertex get_source() const
	{
		return source;
	}
	TVertex get_sink() const
	{
		return sink;
	}
	int get_last(TVertex vertex) const
	{
		return last_edge[vertex];
	}
	int get_last_rev(TVertex vertex) const
	{
		return last_rev_edge[vertex];
	}
	std::vector < int > get_last_edges() const
	{
		return last_edge;
	}
	TFlow get_flow() const
	{
		TFlow current_flow = 0;
		for (int number_edge = last_edge[source]; number_edge != -1; number_edge = prev_edge[number_edge])
		{
			current_flow += edges[number_edge].flow;
		}
		return current_flow;
	}
	friend class IEdgeIterator;
	friend class IOutputEdgeIterator;
	friend class IInputEdgeIterator;
};

typedef int iterator_type;

class IEdgeIterator
{
protected:
	Network* network;
	iterator_type number_edge;
public:
	IEdgeIterator() : network(nullptr), number_edge(-1)
	{
	}
	IEdgeIterator(Network* network_, iterator_type number_edge_) :
		network(network_), number_edge(number_edge_)
	{
	}
	virtual bool valid() const = 0;
	virtual void next() = 0;
	TVertex get_start() const
	{
		return network->edges[number_edge].start;
	}
	int get_number_edge() const
	{
		return number_edge;
	}
	TVertex get_finish() const
	{
		return network->edges[number_edge].finish;
	}
	TFlow get_capacity() const
	{
		return network->edges[number_edge].capacity;
	}
	TFlow get_flow() const
	{
		return network->edges[number_edge].flow;
	}
	NetworkEdge get_edge() const
	{
		return network->edges[number_edge];
	}
	TFlow residual_capacity() const
	{
		return network->edges[number_edge].capacity - network->edges[number_edge].flow;
	}
	bool can_push_flow(TFlow flow) const
	{
		return residual_capacity() - flow >= 0;
	}
	bool saturated() const
	{
		return residual_capacity() == 0;
	}
	void push_flow(TFlow flow)
	{
		network->edges[number_edge].flow += flow;
		network->edges[number_edge ^ 1].flow -= flow;
	}
};

class IOutputEdgeIterator : public IEdgeIterator
{
public:
	IOutputEdgeIterator() : IEdgeIterator()
	{
	}
	IOutputEdgeIterator(Network* network_, TVertex v) : IEdgeIterator(network_, network_->get_last(v))
	{
	}
	bool valid() const
	{
		return network != nullptr && number_edge != -1;
	}
	void next()
	{
		number_edge = network->prev_edge[number_edge];
	}
};

class IInputEdgeIterator : public IEdgeIterator
{
public:
	IInputEdgeIterator() : IEdgeIterator()
	{
	}
	IInputEdgeIterator(Network* network_, TVertex v) : IEdgeIterator(network_, network_->get_last_rev(v))
	{
	}
	bool valid() const
	{
		return network != nullptr && number_edge != -1;
	}
	void next()
	{
		number_edge = network->prev_rev_edge[number_edge];
	}
};

class IAlgorithmFindMaxFlow
{
public:
	virtual TFlow run() = 0;
};

class IAlgorithmPreflow : public IAlgorithmFindMaxFlow
{
protected:
	IAlgorithmPreflow(Network* graph_) : graph(graph_)
	{

	}
	std::vector < TFlow > overflow;
	std::vector < int > height;
	Network* graph;
	static const TFlow INF = 1e9;
	void push(IOutputEdgeIterator& it)
	{
		TVertex u = it.get_start();
		TVertex v = it.get_finish();
		if (overflow[u] > 0 && height[u] == height[v] + 1 && !it.saturated())
		{
			TFlow flow = std::min(it.residual_capacity(), overflow[u]);
			push_flow(flow, it);
			return;
		}
	}
	void relabel(TVertex u)
	{
		if (overflow[u] <= 0)
		{
			return;
		}
		int minheight = INF;
		for (IOutputEdgeIterator it(graph, u); it.valid(); it.next())
		{
			TVertex v = it.get_finish();
			if (!it.saturated())
			{
				if (height[v] < height[u])
				{
					return;
				}
				else
				{
					minheight = std::min(minheight, height[v]);
				}
			}
		}
		height[u] = minheight + 1;
	}
	void push_flow(TFlow flow, IOutputEdgeIterator& it)
	{
		it.push_flow(flow);
		overflow[it.get_start()] -= flow;
		overflow[it.get_finish()] += flow;
	}
};

class AlgorithmPreflow : public IAlgorithmPreflow
{
private:
	TFlow maxflow;
	std::vector < IOutputEdgeIterator > current_edge;
	void discharge(TVertex u)
	{
		while (overflow[u] > 0)
		{
			if (!current_edge[u].valid())
			{
				relabel(u);
				current_edge[u] = { graph, u };
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
public:
	AlgorithmPreflow(Network* graph_) : IAlgorithmPreflow(graph_), maxflow(0)
	{
		overflow.resize(graph->get_vertices_number());
		fill(overflow.begin(), overflow.end(), 0);
		height.resize(graph->get_vertices_number());
		fill(height.begin(), height.end(), 0);
		current_edge.resize(graph->get_vertices_number());
		for (TVertex i = 0; i < graph->get_vertices_number(); i++)
		{
			current_edge[i] = IOutputEdgeIterator(graph, i);
		}
		height[graph->get_source()] = graph->get_vertices_number();
	}
	~AlgorithmPreflow() = default;
	TFlow run()
	{
		TFlow cur_flow = 0;
		for (IOutputEdgeIterator it(graph, graph->get_source()); it.valid(); it.next())
		{
			overflow[graph->get_source()] += it.get_capacity();
			push_flow(it.get_capacity(), it);
		}
		while (true)
		{
			bool can_do_push_or_relabel = false;
			// названия
			for (TVertex i = 0; i < graph->get_vertices_number(); i++)
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
		return graph->get_flow();
	}
};


class IAlgorithmBlockingPath : public IAlgorithmFindMaxFlow
{
protected:
	Network* graph;
	explicit IAlgorithmBlockingPath(Network* graph_) : graph(graph_)
	{
	}
	std::vector<TVertex> find_layers() const
	{
		std::queue<TVertex> q;
		q.push(graph->get_source());
		std::vector<TVertex> level(graph->get_vertices_number());
		level[graph->get_source()] = 1;
		std::vector < NetworkEdge > edges;
		while (!q.empty())
		{
			TVertex u = q.front();
			q.pop();
			for (IOutputEdgeIterator it(graph, u); it.valid(); it.next())
			{
				if (!it.saturated() && (level[it.get_finish()] == 0 || level[it.get_finish()] == level[u] + 1))
				{
					if (level[it.get_finish()] == 0)
					{
						q.push(it.get_finish());
					}
					level[it.get_finish()] = level[u] + 1;
				}
			}
		}
		return level;
	}
	virtual bool find_and_push_blocking_path() = 0;
public:
	TFlow run()
	{
		while (find_and_push_blocking_path())
		{

		}
		return graph->get_flow();
	}
};


class AlgorithmMalhotra : public IAlgorithmBlockingPath
{
private:
	static const TFlow INFFlow = 1e9;
	std::vector < IOutputEdgeIterator > current_edge;
	std::vector < IInputEdgeIterator > current_rev_edge;
	std::vector < TFlow > out_potential, in_potential;
	std::vector < int > bad_edges;
	std::vector < bool > used;
	std::vector < int > layers;
public:
	AlgorithmMalhotra(Network* graph_) : IAlgorithmBlockingPath(graph_)
	{

	}
private:
	TFlow total_potential(TVertex v) const
	{
		if (v == graph->get_source())
		{
			return out_potential[v];
		}
		if (v == graph->get_sink())
		{
			return in_potential[v];
		}
		else
		{
			return std::min(out_potential[v], in_potential[v]);
		}
	}
	template <typename IteratorType>
	void mark_neighbours(TVertex v, std::queue<TVertex>& q, std::vector<TFlow>& potential, bool output)
	{
		for (IteratorType it(graph, v); it.valid(); it.next())
		{
			bad_edges[it.get_number_edge()] = true;
			TVertex number_vertex = (output ? it.get_finish() : it.get_start());
			if (layers[it.get_finish()] == layers[it.get_start()] + 1)
			{
				if (total_potential(number_vertex) != 0)
				{
					potential[number_vertex] -= it.residual_capacity();
					if (total_potential(number_vertex) == 0 && !used[number_vertex])
					{
						used[number_vertex] = true;
						q.push(number_vertex);
					}
				}
			}
		}
	}
	void mark_all_empty_vertices(TVertex v)
	{
		std::queue<TVertex> q;
		q.push(v);
		while (!q.empty())
		{
			TVertex u = q.front();
			q.pop();
			used[u] = true;
			mark_neighbours<IOutputEdgeIterator>(u, q, in_potential, true);
			mark_neighbours<IInputEdgeIterator>(u, q, out_potential, false);
		}
	}
	template < typename IteratorType, typename TFunc, typename TFunc1>
	void push_flow(TVertex u, TFlow flow,  bool direct_push, TFunc edge, TFunc1 potential)
	{
		std::queue< TVertex > q;
		std::vector < TFlow > excess(graph->get_vertices_number());
		std::vector < bool > inqueue(graph->get_vertices_number());
		excess[u] = flow;
		q.push(u);
		inqueue[u] = true;
		while (!q.empty())
		{
			TVertex v = q.front();
			q.pop();
			for (; edge(v).valid(); edge(v).next())
			{
				if (bad_edges[edge(v).get_number_edge()])
				{
					continue;
				}
				TFlow edge_can_flow = edge(v).residual_capacity();
				TVertex u = direct_push ? edge(v).get_finish() : edge(v).get_start();
				if (!excess[u] && !inqueue[u])
				{
					q.push(u);
					inqueue[u] = true;
				}
				TFlow can_flow = std::min(edge_can_flow, excess[v]);
				edge(v).push_flow(can_flow);
				potential(v, u, can_flow);
				excess[u] += can_flow;
				excess[v] -= can_flow;
				if (edge_can_flow > excess[v])
				{
					break;
				}
				else
				{
					bad_edges[edge(v).get_number_edge()] = true;
				}
			}
		}
	}
	void init()
	{
		out_potential.resize(graph->get_vertices_number());
		in_potential.resize(graph->get_vertices_number());
		current_edge.resize(graph->get_vertices_number());
		current_rev_edge.resize(graph->get_vertices_number());
		fill(out_potential.begin(), out_potential.end(), 0);
		fill(in_potential.begin(), in_potential.end(), 0);
		bad_edges.resize(graph->get_edges_number());
		fill(bad_edges.begin(), bad_edges.end(), 0);
		used.resize(graph->get_vertices_number());
		fill(used.begin(), used.end(), false);
	}
	template < typename IteratorType, typename TFunc>
	void make_current_edges(TFunc edge) const
	{
		for (TVertex i = 0; i < graph->get_vertices_number(); i++)
		{
			for (IteratorType it(graph, i); it.valid(); it.next())
			{
				if (layers[it.get_start()] + 1 == layers[it.get_finish()])
				{
					edge(i) = it;
					break;
				}
			}
		}
	}
	void calculate_potential()
	{
		for (TVertex i = 0; i < graph->get_vertices_number(); i++)
		{
			for (IOutputEdgeIterator it(graph, i); it.valid(); it.next())
			{
				if (layers[it.get_start()] + 1 == layers[it.get_finish()])
				{
					in_potential[it.get_finish()] += it.residual_capacity();
					out_potential[it.get_start()] += it.residual_capacity();
				}
				else
				{
					bad_edges[it.get_number_edge()] = true;
				}
			}
		}
	}
	bool find_and_push_blocking_path()
	{
		init();
		layers = find_layers();
		if (layers[graph->get_sink()] == 0)
		{
			return false;
		}
		calculate_potential();
		make_current_edges<IOutputEdgeIterator>([&](TVertex v) -> auto& {return current_edge[v]; });
		make_current_edges<IInputEdgeIterator>([&](TVertex v) -> auto& {return current_rev_edge[v]; });
		for (TVertex i = 0; i < graph->get_vertices_number(); i++)
		{
			if (used[graph->get_source()] || used[graph->get_sink()])
			{
				break;
			}
			bool zero_potential = false;
			TFlow min_potential = INFFlow;
			TVertex number_min;
			for (TVertex v = 0; v < graph->get_vertices_number(); v++)
			{
				if (!used[v] && total_potential(v) == 0)
				{
					mark_all_empty_vertices(v);
					zero_potential = true;
					break;
				}
				if (!used[v])
				{
					if (total_potential(v) < min_potential)
					{
						min_potential = total_potential(v);
						number_min = v;
					}
				}
			}
			if (!zero_potential)
			{
				push_flow<IOutputEdgeIterator>(number_min, min_potential, true, [&](TVertex v) -> auto& {return current_edge[v]; },
					[&](TVertex v, TVertex u, TFlow flow_push) { out_potential[v] -= flow_push; in_potential[u] -= flow_push; });
				push_flow<IInputEdgeIterator>(number_min, min_potential, false, [&](TVertex v) -> auto& {return current_rev_edge[v];  },
					[&](TVertex v, TVertex u, TFlow flow_push) { in_potential[v] -= flow_push; out_potential[u] -= flow_push; });
			}
		}
		return true;
	}
};

struct Data
{
	int n;
	std::vector < int > useless;
	std::vector < std::vector < int > > neighbours;
};

template < typename T, typename FuncT >
std::vector<T> read_vector(std::istream& in, std::size_t vector_size, FuncT func)
{
	std::vector<T> data;
	while (vector_size--)
	{
		T perm;
		in >> perm;
		data.push_back(func(perm));
	}
	return data;
}

template < typename T, typename FuncT>
std::vector<T> read_vector(std::istream& in, FuncT func)
{
	size_t vector_size;
	in >> vector_size;
	return read_vector<T>(in, vector_size, func);
}

template < typename T>
std::vector<T> read_vector(std::istream& in)
{
	return read_vector<T>(in, [](T x){ return x; });
}

Data read(std::istream& in)
{
	Data data;
	data.useless = read_vector<int>(in);
	data.n = data.useless.size();
	data.neighbours.resize(data.n);
	for (int i = 0; i < data.n; i++)
	{
		data.neighbours[i] = read_vector<int>(in, [](int v) {return v - 1; });
	}
	return data;
}

const int INF = 1e6;

void print(std::ostream& out, TFlow ans)
{
	out << ans;
}

TFlow run_algorithm(Data data, int s, int t, std::vector < NetworkEdge> edges, int sum_vertex)
{
	Network* network = new Network(data.n + 2, s, t, edges);
	AlgorithmMalhotra algo(network);
	TFlow max_flow = algo.run();
	delete network;
	return sum_vertex - max_flow;
}

TFlow prepare(Data data)
{
	int s = data.n;
	int t = data.n + 1;
	std::vector < NetworkEdge > edges;
	TFlow sum_vertex = 0;
	for (int i = 0; i < data.n; i++)
	{
		if (data.useless[i] >= 0)
		{
			edges.push_back({ s, i, data.useless[i] });
			sum_vertex += data.useless[i];
		}
		else
		{
			edges.push_back({ i, t, -data.useless[i] });
		}
	}
	for (int v = 0; v < data.n; v++)
	{
		for (int i = 0; i < data.neighbours[v].size(); i++)
		{
			edges.push_back({ v, data.neighbours[v][i], INF });
		}
	}
	return run_algorithm(data, s, t, edges, sum_vertex);
}

void run()
{
	Data data = read(std::cin);
	print(std::cout, prepare(data));
}

int main()
{
#ifdef _DEBUG
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif
	run();
	return 0;
}
