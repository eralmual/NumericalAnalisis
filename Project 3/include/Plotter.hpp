//
// Created by allan on 12/11/18.
//

#ifndef PROYECTO3_PLOTTER_H
#define PROYECTO3_PLOTTER_H

#include <string>
#include <Matrix.hpp>
#include <python2.7/Python.h>
#include <boost/format.hpp>

/**
 * This class calls python to plot the results
 */
class Plotter {
public:
    anpi::Matrix<double> &matrix;
    std::string header, footer;
    int gridScale;
    int height, width;
    bool showFlow;

    Plotter(anpi::Matrix<double> &matrix, int gridScale, int height, int width, bool showFlow)
            : matrix(matrix), gridScale{gridScale},
              height(height), width(width), showFlow(showFlow) {
        header = "import numpy as np\n"
                 "import matplotlib.pyplot as plt\n"
                 "from matplotlib import cm\n"
                 "dpi = 96\n"
                 "margin = 0.05 # (5% of the width/height of the figure...)\n" +
                 (boost::format("xpixels, ypixels = %1%, %2%\n") % width % height).str() +
                 "figsize = (1 + margin) * ypixels / dpi, (1 + margin) * xpixels / dpi\n"
                 "fig = plt.figure(figsize=figsize, dpi=dpi)\n"
                 "ax = fig.add_axes([margin, margin, 1 - 2*margin, 1 - 2*margin])\n";

        footer = "\n"
                 "im=ax.imshow(matrix,interpolation='nearest', cmap=cm.inferno)\n"
                 "plt.colorbar(im, fraction=.046,pad=.04)\n";

    }

    /**
     * Pretty self-explanatory name
     */
    void generatePythonMatrixesAndPlot() {
        Py_Initialize();
        PyRun_SimpleString(header.c_str());
        std::string line("matrix=np.array([\n");

        for (size_t i = 1; i < matrix.rows() - 1; ++i) {
            line += "[";
            for (size_t j = 1; j < matrix.cols() - 1; ++j) {
                line += std::to_string(matrix[i][j]);
                if (j == matrix.cols() - 2) {
                    line += "]";
                } else {
                    line += ",";
                }
            }
            if (i == matrix.rows() - 2) {
                line += "\n])\n";
            } else {
                line += ",\n";
            }
        }

        PyRun_SimpleString(line.c_str());
        PyRun_SimpleString(footer.c_str());
        if (showFlow) {//quiver
            std::cout<<"calculando flujo...\n";
            std::string x = (boost::format("x=np.arange(0,%1%,1)") % (matrix.cols() - 3)).str();
            std::string y = (boost::format("y=np.arange(0,%1%,1)") % (matrix.rows() - 3)).str();
            PyRun_SimpleString(x.c_str());
            PyRun_SimpleString(y.c_str());
            PyRun_SimpleString("xx, yy = np.meshgrid(x, y, sparse=True)");
            PyRun_SimpleString("u=-np.gradient(matrix)[1]");
            PyRun_SimpleString("v=-np.gradient(matrix)[0]");
            PyRun_SimpleString((boost::format(
                    "plt.quiver(xx[::%1%,::%1%],yy[::%1%,::%1%],u[::%1%,::%1%],v[::%1%,::%1%],angles='xy')") %
                                gridScale).str().c_str());
        }
        PyRun_SimpleString("plt.show()");
        Py_Finalize();
    }
};

#endif //PROYECTO3_PLOTTER_H
