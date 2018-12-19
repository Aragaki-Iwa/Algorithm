using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.IO;
using Algorithm;

namespace Algorithm
{
   // class ShortestPath

        //加权有向边的数据类型
        public class DirectedEdge
        {
            private readonly int v;
            private readonly int w;
            private readonly double weight;

            public DirectedEdge(int v, int w, double weight)
            {
                this.v = v; this.w = w; this.weight = weight;
            }

            public double Weight() { return weight; }
            public int From() { return v; }
            public int To() { return w; }
            public string toString() { return string.Format("{0:G}->{1:G} {2:F2}", v, w, weight); }
        }

        //加权有向图 的数据类型
        public class EdgeWeightDigraph
        {
            private readonly int V;//顶点总数
            private int E;
            private NodeList<DirectedEdge>[] adj;

            public EdgeWeightDigraph(int V)
            {
                this.V = V; this.E = 0;
                adj = new NodeList<DirectedEdge>[V];
                for (int i = 0; i < V; i++)
                {
                    adj[i] = new NodeList<DirectedEdge>();
                }
            }
            public EdgeWeightDigraph(string file)
            {
                StreamReader data = new StreamReader(file);
                this.V = int.Parse(data.ReadLine());
                this.E = int.Parse(data.ReadLine());
                adj = new NodeList<DirectedEdge>[V];
                for (int v = 0; v < V; v++)
                {
                    adj[v] = new NodeList<DirectedEdge>();
                }
                string line = data.ReadLine();
                while (line != null)
                {
                    string[] array = line.Split(new char[] { ' ', ',' });
                    int o, d; double w;
                    for (int i = 0; i < array.Length - 1; i++)
                    {
                        o = int.Parse(array[i]);
                        d = int.Parse(array[++i]);
                        w = double.Parse(array[++i]);
                        DirectedEdge e = new DirectedEdge(o, d, w);
                        AddEdge(e);
                    }
                    line = data.ReadLine();
                }

            }

            public int Vnums() { return V; }    public int Enums() { return E; }
            public void AddEdge(DirectedEdge e)
            {
                adj[e.From()].add(e);
                E++;
            }
            public IEnumerable<DirectedEdge> Adj(int v)
            {
                return adj[v];
            }
            public IEnumerable<DirectedEdge> Edges()
            {
                NodeList<DirectedEdge> nodelist = new NodeList<DirectedEdge>();
                for (int v = 0; v < V; v++)
                {
                    foreach (DirectedEdge e in adj[v])
                    {
                        nodelist.add(e);
                    }
                }
                return nodelist;
            }
            public string toString()
            {
                IEnumerable<DirectedEdge> edes = Edges();
                string all = "";
                foreach (DirectedEdge e in edes)
                {
                    all += e.ToString() + "\n";
                }
                return all;
            }
        }

        //基于堆的优先队列
        public class MaxPQ<Key> where Key : IComparable<Key>
        {
            private Key[] heap;
            private int N = 0;

            public MaxPQ(int maxN) { heap = new Key[maxN + 1]; }
            public bool isEmpty() { return N == 0; }
            public int Size() { return N; }
            public void Insert(Key v)
            {
                heap[++N] = v;
                Swim(N);
            }
            public Key DelMax()
            {
                Key max = heap[1];
                Exch(1, N--);
                heap[N + 1] = default(Key);
                Sink(1);
                return max;
            }

            private void Swim(int k)
            {
                while (k > 1 && heap[k].CompareTo(heap[k / 2]) > 0)
                {
                    Exch(k / 2, k);
                    k = k / 2;
                }
            }
            private void Sink(int k)
            {
                while (2 * k <= N)
                {
                    int j = 2 * k;
                    if (j < N && heap[j].CompareTo(heap[j + 1]) < 0) j++;
                    if (heap[k].CompareTo(heap[j]) > 0) break;
                    Exch(k, j);
                    k = j;
                }
            }
            private void Exch(int i, int j)
            { Key r = heap[i]; heap[i] = heap[j]; heap[j] = r; }


        }
        //关联索引的优先队列
        public class IndexMinPQ<Item> where Item : IComparable<Item>
        {
            private Item[] Keys;//元素数组，下标为Index，             值为 元素值
            private int[] heap;//堆，     下标为Index在堆中的位置, 值为Index
            private int[] IndexLocation;      //下标为Index，           值为其在堆中的位置            
            private int N = 0;

            public IndexMinPQ(int maxN)
            {
                Keys = new Item[maxN + 1];
                heap = new int[maxN + 1];
                IndexLocation = new int[maxN + 1];
                for (int i = 0; i < IndexLocation.Length; i++) { IndexLocation[i] = -1; }
            }

            public void Insert(int index, Item item)
            {
                IndexLocation[index] = ++N;
                heap[N] = index;
                Keys[index] = item;

            }
            public Item ItemOf(int index)
            {
                if (Contains(index)) return Keys[index];
                //若无与Index关联的元素
                return default(Item);
            }

            public void Change(int k, Item item)
            {
                Keys[k] = item;
                //改变后的元素大小可能发生变化，或上浮或下沉，为简化处理，均执行
                Swim(IndexLocation[k]); Sink(IndexLocation[k]);
            }
            public bool Contains(int k) { return IndexLocation[k] != -1; }

            public Item Min() { return Keys[heap[1]]; }
            public int MinIndex() { return heap[1]; }
            public int DelMin()
            {
                if (isEmpty()) throw new Exception("队列为空，删不了哦。。。");
                int indexOfMin = heap[1];
                Exch(1, N--);
                Sink(1);
                Keys[indexOfMin] = default(Item);
                IndexLocation[indexOfMin] = -1;
                return indexOfMin;
            }
            public void Delete(int k)
            {
                if (!Contains(k)) { throw new Exception("没有元素与" + k + "关联"); }
                int Location = IndexLocation[k];//此处k代表index，location为它在堆中的位置
                Exch(k, N--);//将k弄到最后，并将堆的长度减1
                // 这里一定要先上浮下沉后再将元素置空，因为swim方法没有N的限制，在没有交换元素的情况下，即删除的就是pq中最后一个元素，如果先置空, 会在Less方法中引发空指针
                // 而sink方法有N的限制，先置空后置空都没有影响，2k <= N会限制它进入循环，避免了空指针
                Swim(Location);
                Sink(Location);
                Keys[k] = default(Item);
                IndexLocation[k] = -1;
            }

            public bool isEmpty() { return N == 0; }
            public int Size() { return N; }

            /***private 是对堆的处理，故其中的k,i,j均代表堆中的位置***/
            /***public 是对元素数组的处理，是基于索引index***/
            //比较的是Keys中元素大小
            private bool Less(int i, int j)
            {
                return Keys[heap[i]].CompareTo(Keys[heap[j]]) < 0;
            }
            private void Exch(int i, int j)
            {
                int index = heap[i];
                heap[i] = heap[j];
                heap[j] = index;
                //还需更新heap的逆数组IndexLocation
                IndexLocation[heap[i]] = j;
                IndexLocation[heap[j]] = i;

            }
            private void Swim(int k)
            {
                while (k > 1 && Less(k, k / 2))
                {
                    Exch(k, k / 2);
                    k = k / 2;//不断上浮
                }
            }
            private void Sink(int k)
            {
                //父结点的最大位置为N/2（因为“二叉”）
                //比较两个子结点，将父结点与较大的互换
                while (2 * k <= N)
                {
                    int j = 2 * k;
                    if (j < N && heap[j].CompareTo(heap[j + 1]) < 0) j++;//确定子结点中大元素 的坐标
                    //判断父结点与所选子结点的大小，若小，停止
                    if (Less(k, j)) break;
                    //否则，互换
                    Exch(k, j);
                    k = j;
                }
            }
        }

        ///最短路
        //Dijkstra
        public class DijkstraSP
        {
            private DirectedEdge[] edgeTo;//中的元素对应的可达顶点构成最短路径树，这棵树即永久标记的点
            private double[] distTo;
            private IndexMinPQ<double> pq;//存放需要被放松的顶点 并确定下一个被放松的顶点；实质就是临时标记的点

            public DijkstraSP(EdgeWeightDigraph G, int s)
            {
                edgeTo = new DirectedEdge[G.Vnums()];
                distTo = new double[G.Vnums()];
                pq = new IndexMinPQ<double>(G.Vnums());
                //将 除起点外的点 的distTo 设为 Infinity
                for (int v = 0; v < G.Vnums(); v++)
                {
                    distTo[v] = double.PositiveInfinity;
                }
                    distTo[s] = 0.0;

                    pq.Insert(s, .0);
                    while (!pq.isEmpty()) Relax(G, pq.DelMin());//删去pq中的min,即min由临时标记变为永久标记
                
            }

            private void Relax(EdgeWeightDigraph G, int v)
            {
                foreach (DirectedEdge e in G.Adj(v))
                {
                    int w = e.To();
                    if (distTo[w] > distTo[v] + e.Weight())
                    {
                        distTo[w] = distTo[v] + e.Weight();
                        edgeTo[w] = e;
                        //若e.to()还未在优先队列中，则Insert it；若在其中 且 优先级需被降低（即找到了比它有更短距离的顶点，也就是先于它被永久标记的顶点），则Change it
                        if (pq.Contains(w)) pq.Change(w, distTo[w]);
                        else pq.Insert(w, distTo[w]);
                    }
                }
            }

            public double DistTo(int v) { return distTo[v]; }
            public bool HasPathTo(int v) { return distTo[v] < double.PositiveInfinity; }
            public IEnumerable<DirectedEdge> PathTo(int v)
            {
                if (!HasPathTo(v)) return null;
                Stack<DirectedEdge> path = new Stack<DirectedEdge>();
                for (DirectedEdge e = edgeTo[v]; e != null; e = edgeTo[e.From()])
                {
                    path.Push(e);
                }
                return path;
            }
        }
        
          //任意顶点对之间的sp
        public class DijkstraAllPairsSP 
        {
            private DijkstraSP[] all;

            public DijkstraAllPairsSP(EdgeWeightDigraph G)
            {
                all = new DijkstraSP[G.Vnums()];
                for (int v = 0; v < G.Vnums(); v++) {
                    all[v] = new DijkstraSP(G, v);
                }
            }
            public IEnumerable<DirectedEdge> s_tPath(int s, int t) {
                return all[s].PathTo(t);
            }
            public double s_tDist(int s, int v) {
                return all[s].DistTo(v);
            }
            
        }
        
        //无环 加权 有向图的 SP Algorithm，时间 与 E+V 正比，与边的 权重是否为负无关
        //先拓扑排序，再松弛
        public class AcycWeDiGSP 
        {
            private DirectedEdge[] edgeTo;
            private double[] distTo;

            public AcycWeDiGSP(EdgeWeightDigraph G, int s) {
                edgeTo = new DirectedEdge[G.Vnums()]; distTo = new double[G.Vnums()];
                for (int v = 0; v < G.Vnums(); v++) {
                    distTo[v] = Double.PositiveInfinity;
                }
                distTo[s] = .0;

                //对 G 先进行拓扑排序
                EWDGraphTopological top = new EWDGraphTopological(G);
                //then, Relax
                foreach (int v in top.Order()) {
                    relax(G, v);
                }
            }
            private void relax(EdgeWeightDigraph G, int v) {
                foreach (DirectedEdge e in G.Adj(v)) {
                    int w = e.To();
                    if (distTo[w] > distTo[v] + e.Weight()) {
                        distTo[w] = distTo[v] + e.Weight();
                        edgeTo[w] = e;
                    }
                }
            }
        }
    /**************/
         //无环 加权 有向图 最长路径
         //①将 所有边的权重 取相反数，则最长==最短路，
         //②将distTo[]的初始值设为 负的最大，并改变relax中 不等式的方向：distTo[w] < distTo[v] + e.weight()
        public class AcycWeDiGLP
        {
            private DirectedEdge[] edgeTo;
            private double[] distTo;

            public AcycWeDiGLP(EdgeWeightDigraph G, int s)
            {
                edgeTo = new DirectedEdge[G.Vnums()]; distTo = new double[G.Vnums()];
                for (int v = 0; v < G.Vnums(); v++)
                {
                    distTo[v] = Double.NegativeInfinity;
                }
                distTo[s] = .0;

                //对 G 先进行拓扑排序
                EWDGraphTopological top = new EWDGraphTopological(G);
                //then, Relax
                foreach (int v in top.Order())
                {
                    relax(G, v);
                }
            }
            private void relax(EdgeWeightDigraph G, int v)
            {
                foreach (DirectedEdge e in G.Adj(v))
                {
                    int w = e.To();
                    if (distTo[w] < distTo[v] + e.Weight())
                    {
                        distTo[w] = distTo[v] + e.Weight();
                        edgeTo[w] = e;
                    }
                }
            }
        }
    //优先级限制的任务调度（多个处理器）
    //每个任务对应3条边：
    //“关键路径”
        public class CPM 
        {
            private int nTask;
            private int S;
            private int T;

            public void CreatNewGraph(string path) {
                StreamReader In = new StreamReader(path);
                int N = int.Parse(In.ReadLine());
                nTask = N;
                S = 2 * N; T = 2 * N;
                EdgeWeightDigraph G = new EdgeWeightDigraph(2 * N + 2);
            }

        }
        
        /// <summary>
        /// 网络流各类问题：最大流，最小费用流，最小费用最大流
        /// </summary>
        public class NetWorkFlow 
        {
            private int V;
            //private int E;
            public int[] prev;
            public double[,] residual;
            public List<ArrayList> AEPath;//增广路集合
            public double maxFlow = 0;

            public double MaxFlow(EdgeWeightDigraph G, int s, int t,int MethodType) //G 中 w 现在代表边的容量
            {
                V = G.Vnums();
                AEPath = new List<ArrayList>();
                //MethodType:1-EK,2-
                if (MethodType == 1) //BFS
                {
                    int u,v;
                    residual = new double[G.Vnums(), G.Vnums()];//残余网络，实际存储的是 两点间的流量
                    prev = new int[V];
                    
                    //初始化
                    for (u = 0; u < G.Vnums(); u++) 
                    {
                        for (v = 0; v < G.Vnums(); v++) 
                        {
                            residual[u, v] = 0;//为0,不存在该边
                        }
                    }
                    foreach (var arc in G.Edges()) 
                    {
                        residual[arc.From(), arc.To()] = arc.Weight();//将其设为边的容量
                    }

                    int p = 0;//记录增广路数量
                    while (BFS(residual, s, t, prev)) //不断寻路，直至找不到
                    {
                        double pathFlow = double.MaxValue;
                        AEPath.Add(new ArrayList());
                       
                        AEPath[p].Add(t);

                        for (v = t; v != s; v = prev[v]) 
                        {
                            u = prev[v];
                            pathFlow = pathFlow < residual[u, v] ? pathFlow : residual[u, v];
                            AEPath[p].Insert(0, u);
                        }
                        
                        AEPath[p].Add(" delta: " + pathFlow);
                        p++;

                        //调整更新
                        for (v = t; v != s; v = prev[v])
                        {
                            u = prev[v];
                            residual[u, v] -= pathFlow;//前向弧流量 +1，即以它的容量 -1 来表示
                            residual[v, u] += pathFlow;
                        }
                        maxFlow += pathFlow;
                    }                    
                }
                return maxFlow;
            }
            private bool BFS(double[,] residual, int s, int t, int[] prev) //寻最短路（没有权重，所以 最短路为s->t 最少的边数）
            {
                bool[] marked = new bool[V];                
                Queue<int> queue = new Queue<int>();

                queue.Enqueue(s);
                marked[s] = true;
                prev[s] = -1;

                while (queue.Count != 0) 
                {
                    int u = queue.Dequeue();

                    for (int v = 0; v < V; v++) 
                    {
                        if (!marked[v] && residual[u, v] > 0) //为访问，且前向弧流量 >0
                        {
                            queue.Enqueue(v);
                            marked[v] = true;
                            prev[v] = u;
                        }
                    }
                }
                //是否访问到了t，即是否存在 s->t 的路
                return marked[t] == true;                
            }
        }

}
