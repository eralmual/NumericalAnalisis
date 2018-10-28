/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * Plotting methods, wrapping Matplotlib (Python)
 * @author: David Badilla S.
 * @email: davbs94@gmail.com
 * @date:   15.07.2017
 */

namespace anpi {

  template <typename T>
  Plot2d<T>::Plot2d(){}

  template <typename T>
  Plot2d<T>::~Plot2d(){}
  
  template <typename T>
  void Plot2d<T>::initialize(int id){
    Py_Initialize();
    PyRun_SimpleString("import matplotlib.pyplot as plt");
    std::string tmp1 = "plt.figure("+std::to_string(id)+")"; 
    PyRun_SimpleString(tmp1.c_str());
    _title = "";
    _xlabel = "";
    _ylabel = "";
    _sizeGrid = 0.f;
  }

  template <typename T>
  void Plot2d<T>::setTitle(const std::string& title){
    _title = title;
    std::string strtitle = "plt.title('"+ _title +"')";
    PyRun_SimpleString(strtitle.c_str());
  }

  template <typename T>
  void Plot2d<T>::setXLabel(const std::string& xlabel){
    _xlabel = xlabel;
    std::string strxlabel = "plt.xlabel('"+ _xlabel +"')";
    PyRun_SimpleString(strxlabel.c_str());
  }
  
  template <typename T>
  void Plot2d<T>::setYLabel(const std::string& ylabel){
    _ylabel = ylabel;
    std::string strylabel = "plt.ylabel('"+ _ylabel +"')";
    PyRun_SimpleString(strylabel.c_str());
  }
  
  template <typename T>
  void Plot2d<T>::setGridSize(const T sizegrid){
    _sizeGrid = sizegrid;
    std::string strgrid = "plt.grid("+std::to_string(_sizeGrid) +")";
    PyRun_SimpleString(strgrid.c_str());
  }
  
  template <typename T>
  void Plot2d<T>::setXRange(const T xi, const T xs){
    std::string strxlim = "plt.xlim("+
      std::to_string(xi) + "," + std::to_string(xs) + ")";
    PyRun_SimpleString(strxlim.c_str());
  }
  
  template <typename T>
  void Plot2d<T>::setYRange(const T yi, const T ys){
    std::string strylim = "plt.ylim("+
      std::to_string(yi) + "," + std::to_string(ys) + ")";
    PyRun_SimpleString(strylim.c_str());
  }
  

  template <typename T>
  void Plot2d<T>::plot(const std::vector<T>& datax,
                       const std::vector<T>& datay,
                       const std::string& label,
                       const std::string& color) {
    std::string xstr = "datax = [";
    std::string ystr = "datay = [";
    std::string pltcmd = "plt.plot(datax,datay";
    if (!label.empty()) {
      pltcmd += ",label='" + label + "'";
    }
    if (!color.empty()) {
      pltcmd += ",color='" + color + "'";
    }
    pltcmd += ")";
    
    for(size_t i = 0; i < datax.size(); i++) {
      if (i == datax.size() - 1){
        xstr.append(std::to_string(datax[i]) + "]");
        ystr.append(std::to_string(datay[i]) + "]");
      } else{
        xstr.append(std::to_string(datax[i]) + ",");
        ystr.append(std::to_string(datay[i]) + ",");
      }
    }

    PyRun_SimpleString(xstr.c_str());
    PyRun_SimpleString(ystr.c_str());
    PyRun_SimpleString(pltcmd.c_str());
    PyRun_SimpleString("plt.legend()");
  }

  template <typename T>
  void  Plot2d<T>::plot(const std::vector<T>& datax,
                        const std::vector<T>& averagey,
                        const std::vector<T>& miny,
                        const std::vector<T>& maxy,
                        const std::string& legend,
                        const std::string& color) {

    // Convert the vectors of data into Python strings
    std::string xstr    = "datax = [";
    std::string avgystr = "avgy = [";
    std::string minystr = "miny = [";
    std::string maxystr = "maxy = [";
    
    char c=',';
    for(size_t i = 0; i < datax.size(); i++) {
      if (i == datax.size()-1) {
        c=']';
      }
      xstr.append(std::to_string(datax[i])   + c);
      avgystr.append(std::to_string(averagey[i])    + c);
      minystr.append(std::to_string(miny[i]) + c);
      maxystr.append(std::to_string(maxy[i]) + c);
    }

    std::string lstr = legend.empty()
      ? ""
      : std::string(",label='") + legend + "'";

    std::string cstr = color.empty()
      ? ""
      : std::string(",color='") + color + "'";

    std::string pltcmd  = "plt.plot(datax,avgy" + lstr + cstr + ",lw=2)";

    std::string fillcmd = "plt.fill_between(datax,miny,maxy" + cstr + ",alpha=0.1)";
    
    // Python lines with the data
    PyRun_SimpleString(xstr.c_str());
    PyRun_SimpleString(avgystr.c_str());
    PyRun_SimpleString(minystr.c_str());
    PyRun_SimpleString(maxystr.c_str());

    // Plot the lines
    PyRun_SimpleString(pltcmd.c_str());
    PyRun_SimpleString(fillcmd.c_str());
    PyRun_SimpleString("plt.legend()");
  }

  
  template <typename T>
  void Plot2d<T>::show(){
    PyRun_SimpleString("plt.show()");
  }
      template<typename T>
      void Plot2d<T>::plot2(const std::vector<T> &datax,
                            const std::vector<T> &datay,
                            const int cols,
                            const int rows,
                            const std::string dir) {
          std::string xstr = "datax = [";
          std::string ystr = "datay = [";
          std::string pltcmd = "plt.plot(datax,datay,color='black')";
          std::string img = "img = plt.imread('";
          img +=  dir;
          img +=  "')";
          std::string fig = "fig, ax = plt.subplots()";
          std::string y = "y = plt.ylim(-0.5,";
          y.append(std::to_string(rows)+"+0.5)");
          std::string x = "x = plt.xlim(-0.5,";
          x.append(std::to_string(cols)+"+0.5)");
          std::string plt = "ax.imshow(img, extent=[0,";
          plt.append(std::to_string(cols)+", 0,"+std::to_string(rows)+"],cmap='gray',interpolation='nearest')");
          for (size_t i = 0; i < datax.size(); i++) {
              if (i == datax.size() - 1) {
                  xstr.append(std::to_string(datax[i]) + "]");
                  ystr.append(std::to_string(datay[i]) + "]");
              } else {
                  xstr.append(std::to_string(datax[i]) + ",");
                  ystr.append(std::to_string(datay[i]) + ",");
              }
          }

          PyRun_SimpleString(xstr.c_str());
          PyRun_SimpleString(ystr.c_str());
          PyRun_SimpleString(img.c_str());
          PyRun_SimpleString(fig.c_str());
          PyRun_SimpleString(x.c_str());
          PyRun_SimpleString(y.c_str());
          PyRun_SimpleString(plt.c_str());
          PyRun_SimpleString(pltcmd.c_str());
          //PyRun_SimpleString("plt.legend()");
      }
      template<typename T>
      void Plot2d<T>::plotQuiver(const anpi::Matrix<float> &U,
                            const anpi::Matrix<float> &V) {
          std::string x_ = "X = [";
          std::string y_ = "Y = [";
          std::string v_ = "V = [";
          std::string u_ = "U = [";
          std::string quiver = "plt.quiver(X, Y, U, V,color='blue', scale_units = 'xy',angles='xy',scale = 1)";
          for (int x = 0; x < U.cols(); x++) {
              for (int y = 0; y < V.rows(); y++) {
                  if (x == U.cols() - 1 && y == V.rows() - 1) {
                      x_ += std::to_string(x) + "]";
                      y_ += std::to_string(-y+V.rows()) + "]";
                      u_ += std::to_string(U[y][x]) + "]";
                      v_ += std::to_string(-V[y][x]) + "]";
                  } else {
                      x_ += std::to_string(x) + ",";
                      y_ += std::to_string(-y+V.rows()) + ",";
                      u_ += std::to_string(U[y][x]) + ",";
                      v_ += std::to_string(-V[y][x]) + ",";
                  }
              }
          }
          PyRun_SimpleString(x_.c_str());
          PyRun_SimpleString(y_.c_str());
          PyRun_SimpleString(u_.c_str());
          PyRun_SimpleString(v_.c_str());
          PyRun_SimpleString(quiver.c_str());
      }
  
} // namespace anpi
