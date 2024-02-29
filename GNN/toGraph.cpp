// g++ --std=c++11 toGraph.cpp gnn_model.cpp -lm -o toGraph.exe
#include <vector>
#include <array>
#include <string>
#include <math.h>
#include "LinkedList.hpp"
#include "PrintHelper.hpp"
#include <iostream>
#include <algorithm>

#include "gnn_model.h"

const unsigned long NMAX=21;
const unsigned long EMAX=42;

// filtering maxes
// float maxAngle = 90; //deg
// int maxSegmLen = 50; // mm
float maxAngle = 15; //deg
float maxSegmLen = 10; // mm
#define USE_FPGA 1
#if (USE_FPGA==1)
float MAX_ANGLE = tan(maxAngle*3.1415/180.);    // 15 in deg
float MAX_SEGM_LEN2 = maxSegmLen* maxSegmLen;   // mm^2 == 10mm;
#endif 

/**
 * @brief 
 * @param offset offset from main array `points`
 * @param nTracks how many tracks to use in one sample
 * @param points 
 * @param nodes[out] shape (NMAX,2)
 * @param truth[out] shape (NMAX,2)
 */
void vector2fixedarray(int offset, int nTracks, std::vector<std::vector<std::vector<float>>>& points, std::vector<std::vector<float>>& nodes, std::vector<std::vector<int>>& truth){
  // Check if maximum is larger than inputted points
  if(points.size() < offset+nTracks)
    offset = 0;

  // initialize nodes and truth
  int nHitID = 0;
  for(int e = offset; e < offset+nTracks; e++){
    // printf("Making from %d -> %d/%d = %d\n",offset, e-offset, nTracks, e);
    int m_id = 0;
    for(int h = 0; h < points.at(e).size(); h++){
      std::vector<float> hit = points[e][h];
      // nodes.push_back({hit[0],hit[2]});
      // truth.push_back({e-offset+1, m_id++});
      nodes[nHitID] = {hit[0],hit[2]}; // can simplify by calling points[e][h][0] and [2] instead
      truth[nHitID] = {e-offset+1, m_id++};
      nHitID++;
      if(nHitID >= NMAX){
        return;
      }
    }
  }
}

/**
 * @brief whether to filter the edge based on a few factors such as distance, angle, and y-distance
 * @param node_i nodes[i]
 * @param node_j nodes[j]
 * @param angle[out] provides the angle between nodes
 * @param dist[out] provides the distance between nodes
 * @returns whether the edge fails filter or not: true if fails, false if passes
 */
bool filterEdge(std::vector<float>& node_i, std::vector<float>& node_j, float& angle, float& dist){

  if(node_i[0] < 0 && node_j[1] < 0)
    return true;
  if(node_i[0] < 0 && node_j[1] < 0)
    return true;
  
  // calculate difference in x,y,z
  std::vector<float> diff;
  for(int q = 0; q < node_i.size(); q++){
    diff.push_back(node_j[q] - node_i[q]);
  }
  // if z position is negative, ignore we only want forward tracks
  if(diff[1] <= 0)
    return true;
  
  // calculate angle and dist between the two nodes
#if (USE_FPGA==0) 
  angle = atan2(diff[0],diff[1]);
  dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]);
  // filter out angle and distance by the maximum values
  if(abs(angle) > maxAngle * (3.1415/180))
    return true;
  if(dist > maxSegmLen)
    return true;
  return false;
#else
  angle = 99.;
  if (diff[1] != 0.) {
    float a = diff[0] / diff[1];
    if (a < 0.)
      a = -a;
    angle = a;
  }
  dist = diff[0] * diff[0] + diff[1] * diff[1]; //hls::sqrt(diff[0]*diff[0]+diff[1]*diff[1]);
  // filter out angle and distance by the maximum values
  if(angle > MAX_ANGLE)
    return true;
  if(dist > MAX_SEGM_LEN2)
    return true;
  return false;
#endif

}


/**
 * @brief Creates edges between nodes
 * @param nodes shape: (events*hits,2) all nodes
 * @param edges[out] shape: (connected_nodes,2) all edges, with connected `senders` and `receivers`
 * @param senders[out] shape: (connected_nodes) all senders, index of node of start of edge
 * @param receivers[out] shape: (connected_nodes) all receivers, index of node of end of edge
 */
void convertToGraph(std::vector<std::vector<float>>& nodes, std::vector<std::vector<float>>& edges, std::vector<int>& senders, std::vector<int>& receivers){

  // iterate over nodes list twice (all combinations)
  int cnt=0;
  for(int i = 0; i < nodes.size(); i++){
    for(int j = 0; j < nodes.size(); j++){
      // if same index (same node) ignore
      if(i == j)
        continue;

      float angle, dist;
      if(filterEdge(nodes[i], nodes[j], angle, dist))
        continue;

      if (cnt>=EMAX) break;

      // populate edges, senders, and receivers
      //edges.push_back({angle,dist});
      //senders.push_back(i);
      //receivers.push_back(j);
      edges.at(cnt)={angle,dist};
      senders.at(cnt)=i;
      receivers.at(cnt)=j;

      std::cout << " cnt= " << cnt << "==>>  i = " << i << " j " << j << " dist = "
		<< dist << " angle = " << angle
	//		<< " edge_counter = " << edge_counter
		<< " tan(MAX_ANGLE) = " << MAX_ANGLE << std::endl;
      cnt++;
    }
  }
}



/**
 * @brief Creates edges between nodes
 * @param nodes shape: (events*hits,2) all nodes
 * @param truth shape: (events*hits,2) all nodes
 * @param edges[out] shape: (connected_nodes,2) all edges
 * @param senders[out] shape: (connected_nodes) all senders, index of node of start of edge
 * @param receivers[out] shape: (connected_nodes) all receivers, index of node of end of edge
 * @param nodes_t[out] shape: (events*hits,1) truthity of nodes, will be all 1's if no ghost hits
 * @param edges_t[out] shape: (connected_nodes,1) truthity of edges, same size as `edges` with 0 or 1, for whether the edge is true
 */
void convertToGraphTest(std::vector<std::vector<float>>& nodes, std::vector<std::vector<int>>& truth, std::vector<std::vector<float>>& edges, std::vector<int>& senders, std::vector<int>& receivers, std::vector<std::vector<float>>& nodes_t, std::vector<std::vector<float>>& edges_t){

  // initialize nodes_t to be all 0's
  nodes_t = std::vector<std::vector<float>>(nodes.size(),std::vector<float>(1,0));

  // iterate over nodes list twice (all combinations)
  for(int i = 0; i < nodes.size(); i++){
    for(int j = 0; j < nodes.size(); j++){
      // if same index (same node) ignore
      if(i == j)
        continue;

      float angle, dist;
      if(filterEdge(nodes[i], nodes[j], angle, dist))
        continue;
          
          
      // populate edges, senders, and receivers
      edges.push_back({angle,dist});
      senders.push_back(i);
      receivers.push_back(j);
      
      // check if the edge is true by the initialized values we made before.
      if(truth[i][0] == truth[j][0] && truth[i][1]+1 == truth[j][1])
        edges_t.push_back({1});
      else
        edges_t.push_back({0});

      nodes_t[i] = {1.};
      nodes_t[j] = {1.};
    }
  }
}

/**
 * 
 * @param tracks[out] array of shape (3,#hits per track)); contains id's of `nodes`
 * @param nodes shape: (events*hits,2) all nodes, shuffled
 * @param edges shape: (connected_nodes,2) all edges
 * @param senders shape: (connected_nodes) all senders, index started of the edge at same index
 * @param receivers shape: (connected_nodes) all receivers, index end of the edge at same index
 * @param nodes_t shape: (events*hits,1) truthity of nodes, will be all 1's if no ghost hits
 * @param edges_t shape: (connected_nodes,1) truthity of edges, same size as `edges` with 0 or 1, for whether the edge is true
 */
void fromGraphToTracks(int nTracks, std::vector<std::vector<int>>& tracks, std::vector<std::vector<float>>& nodes, std::vector<std::vector<float>>& edges, std::vector<int>& senders, std::vector<int>& receivers, std::vector<std::vector<float>>& nodes_t, std::vector<std::vector<float>>& edges_t){

  std::vector<LinkedList*> all_ll;
  bool debug = false;

  for(int eID = 0; eID < edges_t.size(); eID++){
    if(debug) printf("eID: %d\n",eID);
    if(edges_t[eID][0] > 0.5f){
      if(debug) printf(" true\n");
      int sender = senders[eID];
      int receiver = receivers[eID];

      LinkedList* ll = new LinkedList();
      ll->insert_after(sender);
      ll->insert_after(receiver);
      all_ll.push_back(ll);
      if(debug) printf(" appending new: %d->%d\n",sender,receiver);

      for(int i = 0; i < all_ll.size(); i++){
        for(int j = 0; j < all_ll.size(); j++){
          if(i == j) continue;
          if(debug) printf("  checking : %d & %d\n",i,j);
          LinkedList* i_ll = all_ll[i];
          LinkedList* j_ll = all_ll[j];
          if(i_ll->getTail().id == j_ll->getHead().id){
            if(debug) printf("   overlap\n");
            i_ll->resetTail();
            j_ll->resetHead();
            j_ll->remove();
            i_ll->insert_list_after(j_ll);
            all_ll.erase(all_ll.begin()+j);
            j--;
          } else {
            if(debug) printf("   skip\n");
          }
        }
      }

    }
  }

  for(int i = 0; i < all_ll.size(); i++){
    all_ll[i]->resetHead();
    std::vector<int> a;
    for(int j = 0; j < all_ll[i]->size(); j++){
      if(debug) printf("%d:%d | %d\n",i,j,all_ll[i]->getCurrent().id);
      a.push_back(all_ll[i]->getCurrent().id);
      all_ll[i]->next();
    }
    tracks.push_back(a);
  }
  printToPython(tracks);
  printToPython(nodes);
}


#define THRESHOLD 0.5

/**
 * 
 * @param hit_out[out] array of shape (nhits); contains track id of each hit
 * @param nodes shape: (events*hits,2) all nodes, shuffled
 * @param edges shape: (connected_nodes,2) all edges
 * @param senders shape: (connected_nodes) all senders, index started of the edge at same index
 * @param receivers shape: (connected_nodes) all receivers, index end of the edge at same index
 * @param output shape: (connected_nodes) all edges with its connectivity probability
 */
void fromGraph(std::vector<int>& hit_out, int nNodes, std::vector<int>& senders, std::vector<int>& receivers, std::vector<float>& output){

  for(int i = 0; i < hit_out.size(); i++){
    hit_out[i] = 0;
  }

  // printf("sender: %d\n", senders.size());

  int track_id = 1;
  for(int eID = 0; eID < output.size(); eID++){
    if(output[eID] > THRESHOLD){
      int sender = senders[eID];
      int receiver = receivers[eID];

      //printf("s,r: %d,%d\n", sender, receiver);

      int hit_send = hit_out[sender];
      int hit_recv = hit_out[receiver];

      //printf(" s,r: (%d,%d)\n", hit_send, hit_recv);
      //printf(" set to %d\n", track_id);


      if(hit_send == 0 && hit_recv == 0){
	hit_out[sender] = track_id;
	hit_out[receiver] = track_id;
	track_id++;
      } else if(hit_send == 0 && hit_recv != 0){
	hit_out[sender] = hit_recv;
      } else if(hit_send != 0 && hit_recv == 0){
	hit_out[receiver] = hit_send;
      } else { // != 0 && != 0
	for(int hit = 0; hit < nNodes; hit++){
	  if(hit_out[hit] == hit_recv){
	    hit_out[hit] = hit_send;
	  }
	}
      }
      //printf("update:\n");
      //for(int i = 0; i < hit_out.size(); i++){
      //  printf("hit_out[%d] = %d\n", i, hit_out[i]);
      //}
    }
  }
  
  printVector("fromGraph:: tracks: ",hit_out);

}

/**
 * @param x 1D array of all x's
 * @param y 1D array of all y's
 * @param hits_out[out] array of shape (3,#hits per track)); contains id's of `nodes` (`x`,`y` in this case probably)
 */
int doPattern(std::vector<float>& x, std::vector<float>& y, std::vector<int>& hits_out){
  std::vector<std::vector<float>> nodes(NMAX,std::vector<float>(2,-1));
  //std::vector<std::vector<float>> edges;
  //std::vector<int> senders;
  //std::vector<int> receivers;
  std::vector<std::vector<float>> edges(EMAX,std::vector<float>(2,-1));
  std::vector<int> senders(EMAX,0);
  std::vector<int> receivers(EMAX,0);

  printf("x,y size: %d %d\n", x.size(), y.size());
  for(int i = 0; i < std::min(NMAX,x.size()); i++){
    nodes[i][0] = x[i];
    nodes[i][1] = y[i];
  }
  printf("Initialized\n\n");

  printVector("doPattern:: nodes",nodes);
  convertToGraph(nodes, edges, senders, receivers);

  printVector("doPattern:: edges",edges);
  printVector("doPattern:: senders",senders);
  printVector("doPattern:: receivers",receivers);

  printf("====================== gnn_model ===================\n");
  std::vector<float> output;
  gnn_model(nodes, edges, senders, receivers, output);
  printf("====================== end gnn_model ===================\n");

  // std::vector<std::vector<int>> tracks;
  fromGraph(hits_out, nodes.size(), senders, receivers, output);
  return 0;
}
#if 0
int main(){
  int nEvents = 100;
  int nHits = 3;
  /*
  std::vector<std::vector<std::vector<float>>> points = 
  {{{{ 6.59255437, -8.93555313,  5.80234954},
    { 6.6207205,  -8.96406773,  7.30181397},
    { 6.73361779, -9.07844187, 13.31003685},
    { 6.83316992, -9.17938768, 18.60814022},
    { 6.85195454, -9.19843547, 19.60778232},
    { 6.87073516, -9.21748598, 20.60742444},
    { 6.91711939, -9.26456232, 23.07698306}}},
   {{{-1.4460841,   0.70980313,  1.96524883},
    {-1.44776063,  0.71087278,  3.27008381},
    {-1.45234394,  0.71381391,  6.87007969},
    {-1.4527214,   0.71405444,  7.17007936},
    {-1.45792387,  0.71731581, 11.27007476},
    {-1.46212861,  0.71987658, 14.57007109},
    {-1.46644105,  0.7224675,  17.97006736},
    {-1.46999528,  0.7246173,  20.77381182},
    {-1.47228216,  0.72600148, 22.57380983}}},
   {{{ 0.47988469, -3.45631026,  2.56523799},
    { 0.4761946,  -3.4652766,   5.06521919},
    { 0.47487116, -3.46850599,  5.96521242},
    { 0.47147642, -3.47677583,  8.26519505},
    { 0.46566182, -3.49069898, 12.11341656},
    { 0.46000247, -3.50425598, 15.84424039},
    { 0.4541147,  -3.51843827, 19.74421016}}}};
  */
  /*
   std::vector<std::vector<std::vector<float>>> points = 
    { -1.959974, 0,  8.068857}, 
    {-26.438586, 0, 21.114943}, 
    {-25.166266, 0,  8.808744}, 
    {-25.241679, 0,  9.538140}, 
    { -1.831344, 0,  0.065192}, 
    {-26.168207, 0, 18.500206}, 
    { -1.832950, 0,  0.165134}, 
    {-15.225485, 0, 15.711372}, 
    {-14.595705, 0,  4.340377}, 
    {-25.251953, 0,  9.637505}, 
    {-14.397483, 0,  0.760727}, 
    {-24.539594, 0,  2.747473}, 
    { -1.968007, 0,  8.568567}, 
    { -2.051700, 0, 13.765551}, 
    {-14.983134, 0, 11.336239}, 
    {-25.262226, 0,  9.736870}, 
    {-14.391977, 0,  0.661292}, 
    {-14.766441, 0,  7.422853}, 
    { -1.834556, 0,  0.265076}, 
    { -2.120331, 0, 18.023755}, 
    {-26.774119, 0, 24.359648}};
  */

  std::vector<float> x0 = {
     -1.959974,
    -26.438586,
    -25.166266,
    -25.241679,
     -1.831344,
    -26.168207,
     -1.832950,
    -15.225485,
    -14.595705,
    -25.251953,
    -14.397483,
    -24.539594,
     -1.968007,
     -2.051700,
    -14.983134,
    -25.262226,
    -14.391977,
    -14.766441,
     -1.834556,
     -2.120331,
    -26.774119
  };

  std::vector<float> y0 = {
    8.068857, 
    21.114943, 
    8.808744, 
    9.538140, 
    0.065192, 
    18.500206, 
    0.165134, 
    15.711372, 
    4.340377, 
    9.637505, 
    0.760727, 
    2.747473, 
    8.568567, 
    13.765551, 
    11.336239, 
    9.736870, 
    0.661292, 
    7.422853, 
    0.265076, 
    18.023755, 
    24.359648
  };

  std::vector<float> x = {
    6.593,
    6.621,
    6.734,
    6.833,
    6.852,
    6.871,
    6.917,
   -1.446,
   -1.448,
   -1.452,
   -1.453,
   -1.458,
   -1.462,
   -1.466,
   -1.470,
   -1.472,
    0.480,
    0.476,
    0.475,
    0.471,
    0.466
  };

  std::vector<float> y = {
    5.802 ,
    7.302 ,
    13.310,
    18.608,
    19.608,
    20.607,
    23.077,
    1.965 ,
    3.270 ,
    6.870 ,
    7.170 ,
    11.270,
    14.570,
    17.970,
    20.774,
    22.574,
    2.565 ,
    5.065 ,
    5.965 ,
    8.265 ,
    12.113
  };

  std::vector<int> tracks(x.size(), 0);

  doPattern(x, y, tracks);

  // printVector("tracks", tracks);

  for(int i = 0; i < tracks.size(); i++){
    printf(" %d |  %8.2f,%8.2f\n", tracks[i], x[i], y[i]);
  }
  printf("\n\n");


  /*

  printf("\n\n*****************************************************\n");
  printf("***          S T A R T                          *****\n");
  printf("*****************************************************\n\n");

  printVector("points",points);

  std::vector<std::vector<int>> truth(NMAX,std::vector<int>(2,-1));
  std::vector<std::vector<float>> nodes(NMAX,std::vector<float>(2,-1));
  std::vector<std::vector<float>> edges;

  std::vector<std::vector<float>> nodes_t;
  std::vector<std::vector<float>> edges_t;
  std::vector<int> senders;
  std::vector<int> receivers;

  //vector2fixedarray(0, 3, points, nodes, truth);
  printVector("nodes",nodes);
  convertToGraph(nodes, edges, senders, receivers);
  // convertToGraphTest(nodes, truth, edges, senders, receivers, nodes_t, edges_t);

  // Break something in the true data to simulate GNN failing to predict perfect
  // printf("edge_t: %d->%d = %f\n",senders[39], receivers[39], edges_t[39][0]);
  // edges_t[39][0] = 1.;


  printVector("edges",edges);
  printVector("senders",senders);
  printVector("receivers",receivers);

  printf("====================== gnn_model ===================\n");
  gnn_model(nodes, edges, senders, receivers);
  printf("====================== end gnn_model ===================\n");

  std::vector<std::vector<int>> tracks;
  fromGraphToTracks(3, tracks, nodes, edges, senders, receivers, nodes_t, edges_t);
  */
}

#endif
