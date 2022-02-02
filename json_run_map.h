#include "json.hpp"

// Para usar el mapa escribir:
// const Run_map map_of_runs=Run_map("mi_mapa.json");
// Y luego loopear en el objeto json_map como un iterable, por ejemplo:
// for (auto i:map_of_runs.json_map) {
//     if (i["threshold"]>50 && i["threshold"]!=""){cout<<i["run"]<<"  "<< i["threshold"]<<endl;}

class Run_map
{
public:
    nlohmann::json json_map;

    Run_map() {}
    Run_map(string map_file)
    {
        std::ifstream i(map_file);
        i >> json_map;
    }
};
