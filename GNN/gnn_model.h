typedef std::vector<std::vector<float>> f2vec;
typedef std::vector<float> f1vec;

int gnn_model(f2vec nodes,  f2vec edges, std::vector<int> senders, std::vector<int> receivers, f1vec& output);

//#define f2vec(X,Y) hls::vector<hls::vector<t_data,Y>,X>
//#define f1vec(X)   hls::vector<t_data,X>
#define f2vec(X,Y) std::vector<std::vector<float>>
#define f1vec(X)   std::vector<float>
