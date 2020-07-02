//
// Created by Administrator on 3/11/2020.
//

#ifndef GRAND_CHALLENGE_INCREMENTAL_DETECTOR_H
#define GRAND_CHALLENGE_INCREMENTAL_DETECTOR_H

#include <map>
#include <list>
#include <set>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <sys/time.h>

using namespace std;
struct Init_dict
{
    double dbscan_eps;
    int dbscan_min_pts;
    int window_size_n;
    int values_per_second;
    int loss_thresh;
    float temp_eps;
    bool debugging_mode;
    int network_frequency;
};

struct Cluster{
    vector<int> Member_Indices;
    int u;
    int v;
    int cluster_id;
    double Loc;
};

struct Checked_Cluster{
    Cluster c1;
    Cluster c2;
    vector<int> event_interval_t;
};

struct Event
{
    int event_start;
    int event_end;
    int sequenceNumber;
};

struct Feature
{
    double active_power_P;
    double reactive_power_Q;
    int sequenceNumber;
    int badPoint;
};

struct FeatureQueuebuffer
{
    int missingPoint;
    int blocked;
    vector<Feature> buffer;
};

class Incremental_Detector_Disorder {
private:
    double dbscan_eps;
    int dbscan_min_pts;
    int window_size_n;
    int loss_thresh;
    float temp_eps;
    int values_per_second;
    int network_frequency;
    bool debugging_mode;
    bool dbscan_multiprocessing;
    bool is_fitted;
    int MAXIMUM_WINDOW_SIZE;

public:
    vector<Feature> window;
    list<FeatureQueuebuffer> windowBuffer;
    map<int, Feature> mapFeature;
    Cluster cluster0;
    vector<Cluster> mapCluster;
    map<int, Feature> mapFeaturePredict;
    Cluster cluster0Predict;
    vector<Cluster> mapClusterPredict;
    int predictBatchNumber;
    int mapFeatureBlocked;
    int mapFeatureBlockedPos;
    int badPoint;
    int predictCounter;
    Incremental_Detector_Disorder(Init_dict &initDict);
    Feature compute_input_signal(double voltage[1000], double current[1000], int period_length, bool single_sample_mode);
    vector<Event> predict(int batchnumber, Feature X);
    void Incremental_DBSCAN(Feature X);
    void Decremental_DBSCAN(Feature X);
    vector<int> Check_Cluster(Feature X);
    Cluster Check_Cluster0(Feature X);
    double getDis(Feature a, Feature b);
    Checked_Cluster check_event();
    Event backward(Checked_Cluster event_cluster_combination);
    void refresh_window(int MAXIMUM_WINDOW_SIZE);
    void Incremental_DBSCAN_Predict(Feature X);
    vector<int> Check_Cluster_Predict(Feature X);
    Cluster Check_Cluster0_Predict(Feature X);
    void refresh_window_predict(int index);
    int predictEvent(Feature X, int Radius);
};


#endif //GRAND_CHALLENGE_INCREMENTAL_DETECTOR_H
