//
// Created by allan on 11/11/18.
//

#ifndef PROYECTO3_LIEBMANNPARAMS_H
#define PROYECTO3_LIEBMANNPARAMS_H

#include <vector>
#include <IsolationContainer.hpp>
#include <Liebmann.hpp>
#include <Interpolation.hpp>
#include <Matrix.hpp>

/**
 * This class holds the parameters needed to call the anpi::liebmann method
 */
class LiebmannParams {
public:
    int height{}, width{};
    //vector con los perfiles de temperatura
    std::vector<double> topProfile, botProfile, leftProfile, rightProfile;
    anpi::Matrix<double> borders;
    const size_t rowCount=4;

    LiebmannParams() = default;

    /**
     * Checks if a border has no input temperatures, then isolates it
     */
    void verifyIsolation() {
        if(topProfile.empty()){
            container::Isolation::topIsolated=true;
        }
        if(botProfile.empty()){
            container::Isolation::botIsolated=true;
        }
        if(leftProfile.empty()){
            container::Isolation::leftIsolated=true;
        }
        if(rightProfile.empty()){
            container::Isolation::rightIsolated=true;
        }

    }

    /**
     * This method calls Liebmann with the parameters needed
     * @return
     */
    auto calculateLibmann(){
        std::vector<double> topBorderVec,botBorderVec,leftBorderVec,rightBorderVec;

        verifyIsolation();
        auto isolationVector = container::Isolation::convertToIsolationVector();

        borders = anpi::Matrix<double>(rowCount, (size_t)(std::max(height, width)));


        fillBorderVec(topProfile,topBorderVec,(size_t)width);
        fillBorderVec(botProfile,botBorderVec,(size_t)width);

        fillBorderVec(leftProfile,leftBorderVec,(size_t)height);
        fillBorderVec(rightProfile,rightBorderVec,(size_t)height);

        //fill the matrix of frontier conditions
        fillBorderMatrix(topBorderVec, 0);
        fillBorderMatrix(botBorderVec, 1);
        fillBorderMatrix(leftBorderVec, 2);
        fillBorderMatrix(rightBorderVec, 3);

        std::cout<<"calculando liebmann..........\n";
        return anpi::liebmann(borders, (size_t)height, (size_t)width, isolationVector);

    }

private:
    /**
     * This method fills the vector representing each frontier condition depending on the
     * values passed by the user
     * @param profile temperature profile of the user
     * @param frontierVec vector resulting of interpolation according to values given
     * @param size size of the resulting vector
     */
    void fillBorderVec(std::vector<double>& profile,std::vector<double>& frontierVec, const size_t size){
        if(profile.empty()){
            frontierVec =std::vector<double>(size,0);
        }else if(profile.size()==1){
            frontierVec=std::vector<double>(size,profile[0]);
        }
        else if(profile.size()==2){
            //we have two temperatures, so we interpolate linearly
            frontierVec =  anpi::linearInterpolation(0.,double(size-1),profile[0],profile[1],size);

        }else{
            //we have 3 or more temperatures, so we make cubic splines to interpolate
            size_t w= width/(profile.size()-1);
            std::vector<double> xs,ys;

            for (size_t i =0; i< profile.size();++i){
                xs.push_back(i*w);
                ys.push_back(profile[i]);
            }
            frontierVec = anpi::cubicSplinesInterpolation(xs,ys,size);

        }
    }
    /**
     * This method fills the borders of the input matrix for liebmann
     * @param frontier vector containing the temperatures
     * @param direction direction of the vector (1,2,3,4) for (top,bottom,left,right) respectively
     */
    void fillBorderMatrix(std::vector<double> &frontier, size_t direction){
        for (size_t j = 0; j <frontier.size(); ++j){
            borders(direction,j)= frontier[j];
        }
    }


};


#endif //PROYECTO3_LIEBMANNPARAMS_H
