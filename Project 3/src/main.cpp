//
// Created by allan on 06/11/18.
//

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <AnpiConfig.hpp>
#include <Exception.hpp>
#include <boost/algorithm/string.hpp>
#include<boost/algorithm/string/split.hpp>
#include<boost/algorithm/string/compare.hpp>
#include <IsolationContainer.hpp>
#include <Matrix.hpp>
#include <Interpolation.hpp>
#include <Liebmann.hpp>
#include <LiebmannParams.hpp>
#include <Plotter.hpp>


namespace po = boost::program_options;
namespace fs = boost::filesystem;


/**
 * Este método lee el archivo del perfil de temperaturas y retorna un stream
 * de tal forma que si se tiene top = 10 20, la salida sea top=10 top=20,
 * esto para que boost pueda reconocer que se le esta pasando varios valores y los almacene
 * correctamente en un vector
 * @param path ruta del archivo
 */
auto parseConfigFile(const std::string &path) {
    std::cout << "parse config\n";
    std::ifstream file{path};
    std::string line;
    std::string result;
    while (std::getline(file, line)) {

        boost::replace_all(line, "=", " ");
        std::vector<std::string> splitStr;
        boost::algorithm::split(splitStr, line, boost::is_any_of("\t "), boost::token_compress_on);

        std::string parameter = splitStr[0];
        for (auto it = splitStr.begin() + 1; it != splitStr.end(); ++it) {
            result += parameter + "=" + *it + "\n";
        }
    }
    return std::istringstream(result);
}

int main(int argc, const char *argv[]) {
    try {


        LiebmannParams liebmannParams;
        //vector con los perfiles de temperatura
        bool quiet,showFlow;
        int gridScale;

        po::options_description generalOptions("Opciones");
        generalOptions.add_options()
                ("help", "Imprime esta lista de opciones\n")
                ("top,t", po::value<std::vector<double> >(&liebmannParams.topProfile)->multitoken(),
                 "Indica temperatura en borde superior\n")
                ("bottom,b", po::value<std::vector<double> >(&liebmannParams.botProfile)->multitoken(),
                 "Indica temperatura en borde inferior\n")
                ("left,l", po::value<std::vector<double> >(&liebmannParams.leftProfile)->multitoken(),
                 "Indica temperatura en borde izquierdo\n")
                ("right,d", po::value<std::vector<double> >(&liebmannParams.rightProfile)->multitoken(),
                 "Indica temperatura en borde derecho\n")
                ("aislamiento,i", po::value<std::string>(),
                 "<tblr> Aisle los bordes indicados\n \t(t=arriba, b=abajo, l=izquierda, r=derecha)\n \tIndicar uno o más bordes ej '-itblr (todos)'\n")
                ("profile,p", po::value<std::string>(),
                 "Indica el nombre del archivo con perfil termico\n"
                 "Colocarlo en carpeta 'data' de la solución\n")
                ("pixel-horiz,h", po::value<int >(&liebmannParams.width)->default_value(1000),
                 "Número de píxeles horizontales en la solución\n")
                ("pixel-vert,v", po::value<int >(&liebmannParams.height)->default_value(1000),
                 "Número de píxeles verticales en la solución\n")
                ("quiet,q", po::value<bool>(&quiet)->default_value(false),
                 "Desactiva visualizacion\n")
                ("flow,f", po::value<bool>(&showFlow)->default_value(true),
                 "Activa el cálculo de flujo de calor\n")
                ("grid,g", po::value<int>(&gridScale)->default_value(50),
                 "Tamaño de la rejilla de visualización para flujo de calor.\n"
                 "Este valor especifica cuántos pixeles cubre cada celda de la rejilla final.\n"
                 "Por ejemplo si se especifica 5, entonces tanto en x como en y se muestra una flecha del flujo de calor cada 5 pixeles\n");
        //se agrega archivo de perfiles de temperatura
        po::options_description fileOptions{"File"};
        //si se le pasaron los parametros primero, los valores del archivo no son tomados en cuenta
        fileOptions.add_options()
                ("top", po::value<std::vector<double>>(&liebmannParams.topProfile), "top")
                ("bottom", po::value<std::vector<double>>(&liebmannParams.botProfile), "bottom")
                ("left", po::value<std::vector<double>>(&liebmannParams.leftProfile), "left")
                ("right", po::value<std::vector<double>>(&liebmannParams.rightProfile), "right");

        po::variables_map vm;

        store(po::command_line_parser(argc, argv)
                      .options(generalOptions)
                      .extra_parser(container::Isolation::parse_isolation)
                      .run(), vm);

        if (vm["pixel-horiz"].as<int>() < 1) {
            throw anpi::Exception("ancho en pixeles debe ser > 0");
        }
        if (vm["pixel-vert"].as<int>() < 1) {
            throw anpi::Exception("altura en pixeles debe ser > 0");
        }


        if (vm.count("profile")) {
            //se obtiene la ruta completa del archivo deseado
            std::string path = std::string(ANPI_DATA_PATH) + "/" + vm["profile"].as<std::string>();

            bool fileExists = fs::exists(path);
            if (!fileExists) {//se verifica si existe el archivo
                throw anpi::Exception("file not found");
            }

            std::cout << "archivo de temperaturas encontrado: " << path << "\n";
            bool fileIsEmpty = fs::is_empty(path);
            if (fileIsEmpty) {//se verifica si el archivo tiene contenidos
                throw anpi::Exception("file empty");
            }

            auto ifs = parseConfigFile(path);//se formatea el archivo para poder utilizarlo con boost program_options
            if (ifs) {
                po::store(po::parse_config_file(ifs, fileOptions), vm);
            }
        }
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << generalOptions << "\n";//se imprimen las opciones
            return 1;
        }


        auto result = liebmannParams.calculateLibmann();
        std::cout<<"libmann finalizado.\n";

        Plotter plt(result, gridScale, liebmannParams.height, liebmannParams.width, showFlow);

        if(!quiet) {
            std::cout<<"plotendo resultados....\n";
            plt.generatePythonMatrixesAndPlot();
        }
        std::cout<<"end\n";

    }
    catch (const po::error &ex) {
        std::cerr << ex.what() << '\n';
    }
}
