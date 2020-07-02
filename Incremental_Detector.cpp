#include "Incremental_Detector.h"

using namespace std;
vector<int> F1Test;
bool sortFun(Cluster a, Cluster b)
{
    return a.u < b.u;
}

Incremental_Detector::Incremental_Detector( Init_dict &initDict ){
    int VALUES_PER_SECOND = 50;

    Init_dict init_dict;
    init_dict.dbscan_eps = 0.03;
    init_dict.dbscan_min_pts = 2;
    init_dict.window_size_n = 50;
    init_dict.values_per_second = VALUES_PER_SECOND;
    init_dict.loss_thresh = 40;
    init_dict.temp_eps = 0.8;
    init_dict.debugging_mode = false;
    init_dict.network_frequency = 50;

    this->dbscan_eps = init_dict.dbscan_eps;
    this->dbscan_min_pts = init_dict.dbscan_min_pts;
    this->window_size_n = init_dict.window_size_n;
    this->loss_thresh = init_dict.loss_thresh;
    this->temp_eps = init_dict.temp_eps;

    this->values_per_second = init_dict.values_per_second;
    this->network_frequency = init_dict.network_frequency;

    this->debugging_mode = init_dict.debugging_mode;
    this->dbscan_multiprocessing = true;

    this->predictBatchNumber = 0;
    this->predictCounter = 0;

}

Feature Incremental_Detector::compute_input_signal(double voltage[1000], double current[1000], int period_length, bool single_sample_mode){

    Feature X;
    double active_power_P = 0;
    double apparent_power_S = 0;
    double reactive_power_Q = 0;
    double instant_voltage = 0;
    double instant_current = 0;
    double instant_power = 0;
    for(int i = 0; i < period_length; i++){
        instant_voltage += voltage[i]*voltage[i];
        instant_current += current[i]*current[i];
        instant_power += voltage[i] * current[i];
    }
    active_power_P = instant_power/period_length;
    apparent_power_S = sqrt(instant_voltage/period_length)*sqrt(instant_current/period_length);
    reactive_power_Q = sqrt(apparent_power_S * apparent_power_S - active_power_P * active_power_P);

    X.active_power_P = log(active_power_P);
    X.reactive_power_Q =  log(reactive_power_Q);

    return X;
}

Event Incremental_Detector::predict(int batchnumber, Feature X){


    Event event;
    event.event_start = 0;
    event.event_end = 0;
    Checked_Cluster event_cluster_combination;

    bool event_detected = false;
        mapFeature.insert(pair<int, Feature>(X.sequenceNumber, X));
        Incremental_DBSCAN(X);

        if (!cluster0.Member_Indices.empty() && mapCluster.size() >= 2) {
            event_cluster_combination = check_event();
        }
    if(!event_cluster_combination.event_interval_t.empty()){
                Checked_Cluster event_cluster_combination_balanced = event_cluster_combination;
        map<int, Feature>::iterator index;
        for (index=mapFeature.begin(); index!=mapFeature.end(); ++index)
        {
            Feature X_temp = index->second;
            Decremental_DBSCAN(X_temp);
            if (!cluster0.Member_Indices.empty() && mapCluster.size() >= 2) {
                Checked_Cluster event_cluster_combination_balanced_temp = check_event();
                if(!event_cluster_combination_balanced_temp.event_interval_t.empty()){
                    event_cluster_combination_balanced = event_cluster_combination_balanced_temp;
                }else{
                    for(int ind = 0; ind < event_cluster_combination_balanced.event_interval_t.size(); ind++){
                        event_cluster_combination_balanced.event_interval_t[ind] += 1;// the event_interval_t
                    }
                    break;
                }
            }
            else{
                for(int ind = 0; ind < event_cluster_combination_balanced.event_interval_t.size(); ind++){
                    event_cluster_combination_balanced.event_interval_t[ind] += 1;// the event_interval_t
                }
                break;
            }
        }
        int removeFlag = event_cluster_combination.event_interval_t[event_cluster_combination.event_interval_t.size()-1];
        map<int, Feature>::iterator itr = mapFeature.find(removeFlag);
        predictCounter = itr->first;
        mapFeature.erase(mapFeature.begin(),++itr);
        vector<int> tempSetCluster0;
        for(int j = 0 ; j < cluster0.Member_Indices.size(); j++){
            if(cluster0.Member_Indices[j] > removeFlag){
                tempSetCluster0.push_back(cluster0.Member_Indices[j]);
            }
        }
        cluster0.Member_Indices.swap(tempSetCluster0);
        tempSetCluster0.clear();
        for(int i = 0 ; i < mapCluster.size(); i++){
            vector<int> newSet;
            for(int j = 0 ; j < mapCluster[i].Member_Indices.size(); j++){
                if(mapCluster[i].Member_Indices[j] > removeFlag){
                    newSet.push_back(mapCluster[i].Member_Indices[j]);
                }
            }
            mapCluster[i].Member_Indices.swap(newSet);
            newSet.clear();
        }
        vector<Cluster> tempSet;
        for (int i = 0; i < mapCluster.size(); i++) {
            if (!mapCluster[i].Member_Indices.empty()) {
                tempSet.push_back(mapCluster[i]);
            }
        }
        mapCluster.swap(tempSet);
        tempSet.clear();


        event.event_start = event_cluster_combination_balanced.event_interval_t[0];
        event.event_end = event_cluster_combination_balanced.event_interval_t[event_cluster_combination_balanced.event_interval_t.size()-1];
    }
    return event;
}

void Incremental_Detector::Incremental_DBSCAN(Feature X){
    predictBatchNumber = X.sequenceNumber > predictBatchNumber ? X.sequenceNumber : predictBatchNumber;
    if(cluster0.Member_Indices.empty()){
        if(mapCluster.empty()){
            cluster0.Member_Indices.push_back(X.sequenceNumber);
        }
        else{
            vector<int> resultCheck = Check_Cluster(X);
            if(!resultCheck.empty()){
                //put p into clusterN
                mapCluster[resultCheck[0]].Member_Indices.push_back(X.sequenceNumber);
                if(resultCheck.size()>1) {
                    for (int i = 1; i < resultCheck.size(); i++) {
                        //merge
                        mapCluster[resultCheck[0]].Member_Indices.insert(mapCluster[resultCheck[0]].Member_Indices.end(), mapCluster[resultCheck[i]].Member_Indices.begin(),
                                                       mapCluster[resultCheck[i]].Member_Indices.end());
                        mapCluster[resultCheck[i]].Member_Indices.clear();
                    }
                    vector<Cluster> tempSet;
                    for (int i = 0; i < mapCluster.size(); i++) {
                        if (!mapCluster[i].Member_Indices.empty()) {
                            tempSet.push_back(mapCluster[i]);
                        }
                    }
                    mapCluster.swap(tempSet);
                    tempSet.clear();
                }

            }else{
                cluster0.Member_Indices.push_back(X.sequenceNumber);
            }
        }
    }
    else{
        Cluster resultCheckCluster0 = Check_Cluster0(X);
        if(mapCluster.empty()){
            if(resultCheckCluster0.Member_Indices.empty()){
                cluster0.Member_Indices.push_back(X.sequenceNumber);
            }else{
                resultCheckCluster0.Member_Indices.push_back(X.sequenceNumber);
                mapCluster.push_back(resultCheckCluster0);
            }
        }
        else{
            vector<int> resultCheck = Check_Cluster(X);
            if(!resultCheck.empty()){
                mapCluster[resultCheck[0]].Member_Indices.push_back(X.sequenceNumber);

                if(!resultCheckCluster0.Member_Indices.empty()){
                    mapCluster[resultCheck[0]].Member_Indices.insert(mapCluster[resultCheck[0]].Member_Indices.end(),  resultCheckCluster0.Member_Indices.begin(), resultCheckCluster0.Member_Indices.end() );
                }
                if(resultCheck.size() > 1) {
                    for (int i = 1; i < resultCheck.size(); i++) {
                        mapCluster[resultCheck[0]].Member_Indices.insert(mapCluster[resultCheck[0]].Member_Indices.end(), mapCluster[resultCheck[i]].Member_Indices.begin(),
                                                                         mapCluster[resultCheck[i]].Member_Indices.end());
                        mapCluster[resultCheck[i]].Member_Indices.clear();
                    }
                    vector<Cluster> tempSet;
                    for (int i = 0; i < mapCluster.size(); i++) {
                        if (!mapCluster[i].Member_Indices.empty()) {
                            tempSet.push_back(mapCluster[i]);
                        }
                    }
                    mapCluster.swap(tempSet);
                    tempSet.clear();
                }

            }else{
                if(resultCheckCluster0.Member_Indices.empty()){
                    cluster0.Member_Indices.push_back(X.sequenceNumber);
                }else{
                    resultCheckCluster0.Member_Indices.push_back(X.sequenceNumber);
                    mapCluster.push_back(resultCheckCluster0);
                }
            }
        }

    }
}


vector<int> Incremental_Detector::Check_Cluster(Feature X){

    vector<int> resultCheck;
    for(int  i = 0; i < mapCluster.size();i++){
        for(int j =0; j<mapCluster[i].Member_Indices.size();j++){
            Feature X_temp = mapFeature.find(mapCluster[i].Member_Indices[j])->second;
            if(getDis(X, X_temp) <= dbscan_eps) {
                resultCheck.push_back(i);
                break;
            }
        }
    }
    return resultCheck;
}

Cluster Incremental_Detector::Check_Cluster0(Feature X){

    Cluster clusterNew;
    vector<int> newSet;
    for(int i = 0; i<cluster0.Member_Indices.size(); i++){
        Feature X_temp = mapFeature.find(cluster0.Member_Indices[i])->second;
        if(getDis(X, X_temp) <= dbscan_eps) {
            clusterNew.Member_Indices.push_back(X_temp.sequenceNumber);
        }else{
            newSet.push_back(cluster0.Member_Indices[i]);
        }
    }
    cluster0.Member_Indices.swap(newSet);
    newSet.clear();

    return clusterNew;
}

void Incremental_Detector::Decremental_DBSCAN(Feature X){
    int index = X.sequenceNumber;

    if(!cluster0.Member_Indices.empty() && cluster0.Member_Indices[0] == index){
        vector<int> newSet;
        for(int i = 1; i<cluster0.Member_Indices.size(); i++){
                newSet.push_back(cluster0.Member_Indices[i]);
        }
        cluster0.Member_Indices.swap(newSet);
        newSet.clear();
    }else if(!mapCluster.empty()){
        for(int i = 0; i < mapCluster.size(); i++){
            if(mapCluster[i].Member_Indices[0] == index){
                //case 0 cluster size equal to 2
                if(mapCluster[i].Member_Indices.size() == 2){
                    cluster0.Member_Indices.push_back(mapCluster[i].Member_Indices[1]);
                    mapCluster[i].Member_Indices.clear();
                }
                else{
                    vector<int> newSet;
                    mapCluster[i].Member_Indices.swap(newSet);
                    mapCluster[i] = mapCluster.back();
                    mapCluster.pop_back();
                    for(int j = 1; j < newSet.size(); j++){
                        Feature X_temp = mapFeature.find(newSet[j])->second;
                        Incremental_DBSCAN(X_temp);
                    }
                }
                break;
            }
        }
    }

}

void Incremental_Detector::addPoint(Feature X){
}

double Incremental_Detector::getDis(Feature a, Feature b) {

    return sqrt((a.active_power_P-b.active_power_P)*(a.active_power_P-b.active_power_P)+
                (a.reactive_power_Q-b.reactive_power_Q)*(a.reactive_power_Q-b.reactive_power_Q));
}

Checked_Cluster Incremental_Detector::check_event(){
    for(int i = 0; i < mapCluster.size() ; i++){
        sort(mapCluster[i].Member_Indices.begin(), mapCluster[i].Member_Indices.end());
        mapCluster[i].cluster_id = 0;
        mapCluster[i].u = mapCluster[i].Member_Indices[0];
        mapCluster[i].v= mapCluster[i].Member_Indices[mapCluster[i].Member_Indices.size() -1];
        mapCluster[i].Loc = mapCluster[i].Member_Indices.size() / (double) (mapCluster[i].v - mapCluster[i].u + 1);
    }

    vector<Checked_Cluster> checked_clusters_vector;
    vector<int> checked_loc;
    sort(mapCluster.begin(), mapCluster.end(), sortFun);
    for(int i = 0; i < mapCluster.size(); i++){
        if(mapCluster[i].Loc >= 1 - this->temp_eps){
            checked_loc.push_back(i);
        }
    }

    sort(cluster0.Member_Indices.begin(), cluster0.Member_Indices.end());
    for(int i = 0; i < checked_loc.size(); i++){
        for(int j = i+1; j <checked_loc.size(); j++){
            Cluster c1, c2;
            if(mapCluster[checked_loc[i]].u < mapCluster[checked_loc[j]].u){
                c1 = mapCluster[checked_loc[i]];
                c2 = mapCluster[checked_loc[j]];
            }else{
                c2 = mapCluster[checked_loc[i]];
                c1 = mapCluster[checked_loc[j]];
            }

            if(c1.v < c2.u){
                vector<int> event_interval_t;
                Checked_Cluster checkedCluster;
                for(auto index: cluster0.Member_Indices){
                    if(index > c1.v && index < c2.u){
                        event_interval_t.push_back(index);
                    }
                }
                if (!event_interval_t.empty()){
                    checkedCluster.c1 = c1;
                    checkedCluster.c2 = c2;
                    checkedCluster.event_interval_t = event_interval_t;
                    checked_clusters_vector.push_back(checkedCluster);
                }
            }
        }
    }

    Checked_Cluster checkedCluster;
    if(!checked_clusters_vector.empty()) {
        int min_index = 0;
        double event_model_loss = -1;
        for (int i = 0; i < checked_clusters_vector.size(); i++) {
        int lower_event_bound_u = checked_clusters_vector[i].event_interval_t[0] - 1;
        int upper_event_bound_v = checked_clusters_vector[i].event_interval_t[checked_clusters_vector[i].event_interval_t.size()-1] - 1;

            int event_model_loss_temp = 0;
            for (int temp: checked_clusters_vector[i].c1.Member_Indices) {
                if (temp > lower_event_bound_u) {
                    event_model_loss_temp++;
                }
            }
            for (int temp: checked_clusters_vector[i].c2.Member_Indices) {
                if (temp < upper_event_bound_v) {
                    event_model_loss_temp++;
                }
            }
            if (event_model_loss > event_model_loss_temp || event_model_loss == -1) {
                event_model_loss = event_model_loss_temp;
                min_index = i;
            }
        }
        if (event_model_loss <= this->loss_thresh) {
            checkedCluster = checked_clusters_vector[min_index];
            checkedCluster.c1.cluster_id = 1;
        }
    }
    return checkedCluster;
}

void Incremental_Detector::refresh_window(int MAXIMUM_WINDOW_SIZE){
    if(mapFeature.size() > MAXIMUM_WINDOW_SIZE){
        predictCounter = mapFeature.rbegin()->first;
        mapFeature.clear();
        cluster0.Member_Indices.clear();
        mapCluster.clear();
    }
}