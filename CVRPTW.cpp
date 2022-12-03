#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <climits>
#include <queue>
thread_local std::mt19937 gen{ std::random_device{}() };

template<typename T>
T random(T min, T max) {
    return std::uniform_int_distribution<T>{min, max}(gen);
}

using namespace std;
using namespace std::chrono;

struct Customer {
    int ID, X, Y, Demand, ReadyTime, DueTime, Service;
    Customer();
    ~Customer() = default;
};
Customer::Customer() {
    ID = 0;
    X = 0;
    Y = 0;
    Demand = 0;
    ReadyTime = 0;
    DueTime = 0;
    Service = 0;
}
struct Track {
    int Capacity;
    vector<Customer> Route;
    Track(int TrackCapacity);
    ~Track() = default;
};
Track::Track(int TrackCapacity) {
    Capacity = TrackCapacity;
}
struct Result {
    int Track = 0;
    long double Time = 0;
    Result(int x, long double y);
};
Result::Result(int x, long double y) {
    Track = x;
    Time = y;
}

bool ReadFile(string& Filename, vector<Customer>& CustomersVector, int& TrackCapacity) {
    string Tmp[15];
    string Check[15] = { "VEHICLE", "NUMBER", "CAPACITY", "CUSTOMER", "CUST", "NO.", "XCOORD.", "YCOORD.", "DEMAND",
                        "READY", "TIME", "DUE", "DATE", "SERVICE", "TIME" };
    string ProblemID;
    int NumberOfTrack;

    ifstream File;
    File.open(Filename);
    if (!File.good()) {
        cout << "File open error" << endl;
        return false;
    }
    else {
        File >> ProblemID;

        for (int i = 0; i < 3; i++) {
            File >> Tmp[i];
            if (Tmp[i] != Check[i]) {
                cout << "Invalid data format" << endl;
                return false;
            }
        }
        File >> NumberOfTrack >> TrackCapacity;
        for (int i = 3; i < 15; i++) {
            File >> Tmp[i];
            if (Tmp[i] != Check[i]) {
                cout << "Invalid data format" << endl;
                return false;
            }
        }

        while (true) {
            Customer ToAdd;
            File >> ToAdd.ID >> ToAdd.X >> ToAdd.Y >> ToAdd.Demand >> ToAdd.ReadyTime >> ToAdd.DueTime
                >> ToAdd.Service;
            if (File.eof()) break;
            CustomersVector.push_back(ToAdd);
        }
    }
    File.close();
    return true;
};
bool ReadingFromFileInParts(string& Filename, vector<Customer>& CustomersVector, int& TrackCapacity, int CustomersNumber) {
    string Tmp[15];
    string Check[15] = { "VEHICLE", "NUMBER", "CAPACITY", "CUSTOMER", "CUST", "NO.", "XCOORD.", "YCOORD.", "DEMAND",
                        "READY", "TIME", "DUE", "DATE", "SERVICE", "TIME" };
    string ProblemID;
    int NumberOfTrack;

    ifstream File;
    File.open(Filename);
    if (!File.good()) {
        cout << "File open error" << endl;
        return false;
    }
    else {
        File >> ProblemID;

        for (int i = 0; i < 3; i++) {
            File >> Tmp[i];
            if (Tmp[i] != Check[i]) {
                cout << "Invalid data formath" << endl;
                return false;
            }
        }
        File >> NumberOfTrack >> TrackCapacity;
        for (int i = 3; i < 15; i++) {
            File >> Tmp[i];
            if (Tmp[i] != Check[i]) {
                cout << "Invalid data format" << endl;
                return false;
            }
        }

        for (int i = 0; i <= CustomersNumber; i++) {
            Customer ToAdd;
            File >> ToAdd.ID >> ToAdd.X >> ToAdd.Y >> ToAdd.Demand >> ToAdd.ReadyTime >> ToAdd.DueTime
                >> ToAdd.Service;
            if (File.eof()) break;
            CustomersVector.push_back(ToAdd);
        }
    }
    File.close();
    return true;
};
bool TestCondition(vector<Customer>& CustomersVector, int& TrackCapacity, long double** Matrix) {
    long double value;
    for (int i = 1; i < CustomersVector.size(); i++) {
        value = 0;
        if (CustomersVector[i].Demand > TrackCapacity) return false;                                                  //czy zapotrzebowanie klienta nie jest większe od ładowności ciężarówki?
        if (Matrix[0][i] > CustomersVector[i].DueTime) return false;                                                    //czy do każdego sklepu można dojechać przed jego zamknieciem?
        value += max(Matrix[0][i], (long double)CustomersVector[i].ReadyTime);
        value += CustomersVector[i].Service;
        value += Matrix[i][0];
        if (value > CustomersVector[0].DueTime) return false;                                                           //czy zawsze zdaży się spowrotem do magazynu?
    }
    return true;
}

long double DistanceEvaluator(Customer& A, Customer& B) {
    long double value = sqrt(
        (long double)pow((A.X - B.X), 2) + (long double)pow((A.Y - B.Y), 2));
    return value;
}
void DistanceMatrix(long double** Matrix, vector<Customer>& CustomersVector) {
    for (int i = 0; i < CustomersVector.size(); i++)
        for (int j = 0; j < CustomersVector.size(); j++) {
            if (i == j) Matrix[i][j] = 0;
            else Matrix[i][j] = DistanceEvaluator(CustomersVector[i], CustomersVector[j]);
        }
    return;
}
long double ServiceBeginning(vector<Customer>& Route, int No) {
    if (No == 0) return 0;
    long double value =
        ServiceBeginning(Route, No - 1) + Route[No - 1].Service + DistanceEvaluator(Route[No - 1], Route[No]);
    if (value > Route[No].ReadyTime) return value;
    else return Route[No].ReadyTime;
}                                                     //moment rozpoczęcia obsługi danego punktu
long double SolutionValue(vector<Track>& TrackVector) {
    long double value = 0;
    for (int i = 0; i < TrackVector.size(); i++) {
        int LastCustomer = TrackVector[i].Route.size() - 1;
        value += ServiceBeginning(TrackVector[i].Route, LastCustomer) +
            TrackVector[i].Route[LastCustomer].Service +
            DistanceEvaluator(TrackVector[i].Route[LastCustomer], TrackVector[i].Route[0]);
    }
    return value;
}
void SolutionPrinting(vector<Track>& TrackVector) {
    cout << TrackVector.size() << " " << SolutionValue(TrackVector) << endl;
}
bool SaveFile(vector<Track>& TrackVector, string output_filename) {
    ofstream File;
    File.open(output_filename);
    if (!File.good()) {
        cout << "Error opening file for writing" << endl;
        return false;
    }
    else {
        File << TrackVector.size() << " " << setprecision(16) << SolutionValue(TrackVector) << "\n";
        for (int i = 0; i < TrackVector.size(); i++) {
            for (int j = 1; j < TrackVector[i].Route.size(); j++) {
                File << TrackVector[i].Route[j].ID << " ";
            }
            File << "\n";
        }
    }
    File.close();
    return true;
}

void GreedySearch(vector<Customer>& CustomersVector, vector<Track>& TrackVector, int& TrackCapacity, long double** Matrix) {
    bool* VisitedCustomers = new bool[CustomersVector.size()];
    for (int i = 0; i < CustomersVector.size(); i++) VisitedCustomers[i] = false;

    bool Done = false;
    int CurrentVehicle = 0;
    int CustomerNo = 0;

    while (!Done) {
        long double AvailableTime = CustomersVector[0].DueTime;

        int AvailableCapacity = TrackCapacity;
        long double CurrentTime = 0;
        int CurrentPosition = 0;

        TrackVector.push_back(Track(TrackCapacity));
        TrackVector[CurrentVehicle].Route.push_back(CustomersVector[0]);

        bool Possible = true;
        long double BestValue;
        long double PossibleValue = 0;
        int CustomerToVisit = 0;
        while (Possible) {                                                                                              //szukanie najblizszego wierzcholka
            BestValue = INT_MAX;
            Possible = false;
            for (int i = 1; i < CustomersVector.size(); i++) {
                if (VisitedCustomers[i]) continue;                                                                      //czy klient był już odwiedzony?
                if (CustomersVector[i].Demand > AvailableCapacity) continue;                                            //czy ciężarówka ma dość ładunku, by obsłużyć danego klienta?
                if (CurrentTime + Matrix[CurrentPosition][i] > CustomersVector[i].DueTime) continue;                    //czy ciężarówka zdąży przyjechać przed zamknięciem sklepu?
                if (max(CurrentTime + Matrix[CurrentPosition][i], (long double)CustomersVector[i].ReadyTime) +
                    CustomersVector[i].Service + Matrix[i][0] > AvailableTime)
                    continue;                                                                                           //czy ciężarówka zdąży wrócić do magazynu przed jego zamknęciem?

                Possible = true;
                PossibleValue = max(Matrix[CurrentPosition][i],
                    (long double)CustomersVector[i].ReadyTime - CurrentTime);
                if (PossibleValue < BestValue) {                                                                        //warunek sprawdza czy warto jechac do i-tego wierzcholka
                    BestValue = PossibleValue;
                    CustomerToVisit = i;
                }
            }
            if (!Possible) break;

            TrackVector[CurrentVehicle].Route.push_back(CustomersVector[CustomerToVisit]);
            VisitedCustomers[CustomerToVisit] = true;
            CustomerNo++;
            CurrentTime += max(Matrix[CurrentPosition][CustomerToVisit],
                (long double)CustomersVector[CustomerToVisit].ReadyTime - CurrentTime) +
                CustomersVector[CustomerToVisit].Service;
            AvailableCapacity -= CustomersVector[CustomerToVisit].Demand;                                                   //zapas pozostałego ładunku w ciężarówce
            CurrentPosition = CustomerToVisit;
        }

        CurrentVehicle++;
        Done = true;                                                                                                    //warunek zakonczenia petli po odwiedzeniu wszystkich wierzcholkow
        for (int i = 1; i < CustomersVector.size(); i++) {
            if (!VisitedCustomers[i])
                Done = false;
        }
    }

    delete[] VisitedCustomers;
    return;
}




bool fitness(vector<Track>& Base, vector<Track>& Candidate) {                    //funkcja zwraca True, jesli Candidate jest lepszy od Base
    if (Candidate.size() < Base.size())
        return true;
    else if (Candidate.size() == Base.size() && SolutionValue(Candidate) < SolutionValue(Base))
        return true;
    else
        return false;
}

bool member(queue< vector<Track> >& Queue, vector<Track>& Member) {              //funkja zwraca True, jesli Member nalezy do kolejki Queue

    vector <Track> temp;
    int queue_size = Queue.size();
    bool result = false;
    bool conditionA, conditionB;
    for (int i = 0; i < queue_size; i++)
    {
        temp.clear();
        temp = Queue.front();
        Queue.pop();
        conditionA = true, conditionB = true;
        if (result == false) {
            if (temp.size() != Member.size())
            {
                Queue.push(temp); continue;
            }

            for (int j = 0; j < Member.size(); j++) {
                if (Member[j].Route.size() != temp[j].Route.size() || conditionB == false)
                {
                    conditionA = false; break;
                }
                for (int k = 0; k < Member[j].Route.size(); k++) {
                    if (Member[j].Route[k].ID != temp[j].Route[k].ID)
                    {
                        conditionB = false; break;
                    }
                }
            }
            if (conditionA && conditionB)
                result = true;
        }
        Queue.push(temp);
    }
    return result;
}

int track_total_capacity(vector<Customer>& CustomersVector, Track& truck)           // ile ciezarowka ladunku zuzywa
{
    int capacity = 0;
    for (int i = 1; i < truck.Route.size(); i++)
        capacity += truck.Route[i].Demand;
    return capacity;
}

bool validate_route(vector<Customer>& CustomersVector, Track& truck, long double** Matrix)   //funkcja zwraca True jak droga jest poprawna
{
    int CurrentPosition, PreviousPosition = 0;
    long double CurrentTime = 0;
    for (int i = 1; i < truck.Route.size(); i++)
    {
        CurrentPosition = truck.Route[i].ID;
        CurrentTime += max(Matrix[PreviousPosition][CurrentPosition],
            (long double)CustomersVector[CurrentPosition].ReadyTime - CurrentTime) +
            CustomersVector[CurrentPosition].Service;
        if (CurrentTime - CustomersVector[CurrentPosition].Service > CustomersVector[CurrentPosition].DueTime)
            return false;
        PreviousPosition = CurrentPosition;
    }
    if (CurrentTime + Matrix[CurrentPosition][0] > CustomersVector[0].DueTime)
        return false;
    else
        return true;
}

bool Permutation(vector<Customer>& CustomersVector, vector<Track>& Primary, long double** Matrix)            //zwraca True jesli rozwiazanie sie zmienilo
{
    Track TempTruck(Primary[0].Capacity);

    int vehicle_n0 = random(0, (int)(Primary.size() - 1));
    int vehicle_route_size = Primary[vehicle_n0].Route.size();
    int insert_node_n0 = random(1, vehicle_route_size - 1);

    int insert_node_CID = Primary[vehicle_n0].Route[insert_node_n0].ID;
    for (int i = 0; i < Primary.size(); i++)
    {
        if (i == vehicle_n0)                    //aby nie wlozyc wierzcholka w to samo miejsce
            continue;
        if (track_total_capacity(CustomersVector, Primary[i]) + CustomersVector[insert_node_CID].Demand > Primary[i].Capacity) //czy starczy miejsca w ciezarowce
            continue;

        TempTruck.Route.clear();
        TempTruck.Route.push_back(CustomersVector[0]);
        bool pushed = false;
        long double PossibleTime, CurrentTime = 0;
        int CurrentPosition, PreviousPosition = 0;
        for(int j = 1; j < Primary[i].Route.size(); j++)
        {
            CurrentPosition = Primary[i].Route[j].ID;

            PossibleTime = CurrentTime + max(Matrix[PreviousPosition][insert_node_CID],
                (long double)CustomersVector[insert_node_CID].ReadyTime - CurrentTime) +
                CustomersVector[insert_node_CID].Service;
            CurrentTime += max(Matrix[PreviousPosition][CurrentPosition],
                (long double)CustomersVector[CurrentPosition].ReadyTime - CurrentTime) +
                CustomersVector[CurrentPosition].Service;
            if (PossibleTime < CurrentTime && !pushed)
            {
                TempTruck.Route.push_back(CustomersVector[insert_node_CID]);
                pushed = true;
                CurrentTime = PossibleTime;
                CurrentTime += max(Matrix[insert_node_CID][CurrentPosition],
                    (long double)CustomersVector[CurrentPosition].ReadyTime - CurrentTime) +
                    CustomersVector[CurrentPosition].Service;
            }

            TempTruck.Route.push_back(Primary[i].Route[j]);

            PreviousPosition = CurrentPosition;
        }
        if (pushed && validate_route(CustomersVector, TempTruck, Matrix))
        {

            Primary[i].Route.clear();
            for (int a = 0; a < TempTruck.Route.size(); a++)
                Primary[i].Route.push_back(TempTruck.Route[a]);

            if (vehicle_route_size == 2)                    //jesli ciezarowka jechala tylko do jednego klienta to trzeba ja usunac
                Primary.erase(Primary.begin() + vehicle_n0);
            else                                            // a jak nie to trzeba wyrzucic wierzcholek z jej trasy
                Primary[vehicle_n0].Route.erase(Primary[vehicle_n0].Route.begin() + insert_node_n0);

            return true;
        }

    }
    return false;
}


void searchNeighbors(vector<Customer>& CustomersVector, vector< vector<Track> >& Neighborhood, vector<Track>& Candidate, long double** Matrix, int Neighbors, int TabuMoves)   //szuka sasiedztwa rozwiazania
{
    vector<Track> temp_vector;

    for (int a = 0; a < Neighbors; a++)
    {
        temp_vector.clear();
        for (int i = 0; i < Candidate.size(); i++)    //przepisanie Kandidata do tempa
            temp_vector.push_back(Candidate[i]);

        for (int b = 0; b < TabuMoves; b++)
        {
           
            Permutation(CustomersVector, temp_vector, Matrix);
            
        }
        Neighborhood.push_back(temp_vector);
    }
    return;
}

void Neighbors_sorting(vector< vector<Track> >& SNeighborhood)
{
    int index, min_vehicles;
    long double min_result;
    vector< vector<Track> > SNeighborhood_temp;
    for (int i = 0; i < SNeighborhood.size(); i++)
        SNeighborhood_temp.push_back(SNeighborhood[i]);

    int iterations = SNeighborhood_temp.size();
    SNeighborhood.clear();
    for (int j = 0; j < iterations; j++)
    {
        min_result = 1000000;
        min_vehicles = 10000000;
        for (int i = 0; i < SNeighborhood_temp.size(); i++)
        {
            if (SNeighborhood_temp[i].size() < min_vehicles)
            {
                min_vehicles = SNeighborhood_temp[i].size();
                min_result = SolutionValue(SNeighborhood_temp[i]);
                index = i;
            }
            else if (SNeighborhood_temp[i].size() == min_vehicles && SolutionValue(SNeighborhood_temp[i]) < min_result)
            {
                min_result = SolutionValue(SNeighborhood_temp[i]);
                index = i;
            }
        }
        SNeighborhood.push_back(SNeighborhood_temp[index]);
        SNeighborhood_temp.erase(SNeighborhood_temp.begin() + index);
    }

    SNeighborhood_temp.clear();
    return;
}

void Tabu(vector<Customer>& CustomersVector, vector<Track>& TrackVector, int& TrackCapacity, long double** Matrix, double& TimeLimit, int TabuSize, int Neighbors, int TabuMoves) {
    high_resolution_clock::time_point t1, t2;
    duration<double> Timer = (duration<double>) 0;
    t1 = high_resolution_clock::now();

    queue< vector<Track> > TabuList;
    vector< vector<Track> > SNeighborhood;

    vector<Track> SBest;
    for (int i = 0; i < TrackVector.size(); i++) {
        SBest.push_back(TrackVector[i]);
    }
    vector<Track> BestCandidate;
    for (int i = 0; i < TrackVector.size(); i++) {
        BestCandidate.push_back(TrackVector[i]);
    }

    TrackVector.clear();
    TabuList.push(TrackVector);

    bool Clock = false;
    while (!Clock)
    {
        t2 = high_resolution_clock::now();
        Timer = duration_cast<duration<double>>(t2 - t1);
        if (Timer.count() > TimeLimit) {
            Clock = true;
        }

        SNeighborhood.clear();
        searchNeighbors(CustomersVector, SNeighborhood, BestCandidate, Matrix, Neighbors, TabuMoves);
        Neighbors_sorting(SNeighborhood);

        for (int i = 0; i < SNeighborhood.size(); i++) {
            if (!member(TabuList, SNeighborhood[i])) {
                BestCandidate.clear();
                for (int j = 0; j < SNeighborhood[i].size(); j++) {
                    BestCandidate.push_back(SNeighborhood[i][j]);
                }
                break;
            }
        }

        if (fitness(SBest, BestCandidate)) {
            SBest.clear();
            for (int j = 0; j < BestCandidate.size(); j++) {
                SBest.push_back(BestCandidate[j]);
            }
        }

        TabuList.push((BestCandidate));

        if (TabuList.size() > TabuSize)
            TabuList.pop();
    }

    for (int i = 0; i < SBest.size(); i++) {
        TrackVector.push_back(SBest[i]);
    }

    return;
}




Result TabuStart(vector<Customer>& CustomersVector, vector<Track>& BestVector, int& TrackCapacity, long double** DistancesMatrix, double& TimeLimit, int TabuSize, int Neighbors, int TabuMoves) {

    long double TimeResult, BestTime;
    int TracksResult, BestVehicles;

    GreedySearch(CustomersVector, BestVector, TrackCapacity, DistancesMatrix);
    
    BestTime = SolutionValue(BestVector);
    BestVehicles = (int)BestVector.size();
    Tabu(CustomersVector, BestVector, TrackCapacity, DistancesMatrix, TimeLimit, TabuSize, Neighbors, TabuMoves);
    TimeResult = SolutionValue(BestVector);
    TracksResult = (int)BestVector.size();

    Result value = Result(BestVehicles, BestTime);
    return value;
}

void Initialization(string input, string output, double time, int TabuSize, int Neighbors, int TabuMoves) {
    double TimeLimit = time;
    string Filename;
    vector<Customer> CustomersVector;
    vector<Track> TrackVector;
    int TrackCapacity;

    Filename = input;
    if (!ReadFile(Filename, CustomersVector, TrackCapacity))exit(-1);
    
    long double** DistancesMatrix = new long double* [CustomersVector.size()];
    for (int i = 0; i < CustomersVector.size(); i++) DistancesMatrix[i] = new long double[CustomersVector.size()];
    DistanceMatrix(DistancesMatrix, CustomersVector);

    if (!TestCondition(CustomersVector, TrackCapacity, DistancesMatrix)) {
        cout << "No solution" << endl;
        ofstream File;
        File.open(output);
        if (!File.good()) {
            cout << "File reading error" << endl;
            exit(-1);
        }
        else {
            File << -1;
        }
    }
    else {
        TabuStart(CustomersVector, TrackVector, TrackCapacity, DistancesMatrix, TimeLimit, TabuSize, Neighbors, TabuMoves);
        SolutionPrinting(TrackVector);
        SaveFile(TrackVector, output);
    }

    for (int i = 0; i < CustomersVector.size(); i++)
        delete[] DistancesMatrix[i];
    delete[] DistancesMatrix;
}

int main(int argc, char* argv[]) {
    double TimeLimit = 60;
    int TabuSize = 15;
    int Neighbors = 60;
    int TabuMoves = 3;
    string input_filename = argv[1];
    string output_filename = argv[2];
     if (argc >= 4)
        TimeLimit = atof(argv[3]);
    Initialization(input_filename, output_filename, TimeLimit, TabuSize, Neighbors, TabuMoves);

    return 0;
}