#include<iostream>
#include<math.h>
#include<cmath>
#include<vector>
#include<fstream>
#include<algorithm>
#include<sstream>
#define EarthRadius 6371.0
#define epsilon 0.00001
#define deg2rad M_PI/180
#define rad2deg 180/M_PI

class Rays
{
      public:
	    double lat0;
	    double lon0;
	    double lat1;
	    double lon1;

      public:
	    void print(std::ofstream & out)
	    {
		  out << lon0 << "\t";
		  out << lat0 << "\t";
		  out << lon1 << "\t";
		  out << lat1 << "\t";
	    }
	    double length()
	    {
		  double lon0r = deg2rad*lon0;
      		  double lon1r = deg2rad*lon1;
      		  double lat0r = deg2rad*lat0;
      		  double lat1r = deg2rad*lat1;
      		  double dlat = lat1r - lat0r;
      		  double dlon = lon1r - lon0r;
		  double a = sin(dlat/2.0) * sin(dlat/2.0) + cos(lat0r) * cos(lat1r) * sin(dlon/2.0) * sin(dlon/2.0);
		  double c = 2.0 * atan2(sqrt(a),sqrt(1-a));
		  return EarthRadius * c;
	    }
	    double initialAzimuth()
	    {
		  double lon0r = deg2rad*lon0;
      		  double lon1r = deg2rad*lon1;
      		  double lat0r = deg2rad*lat0;
      		  double lat1r = deg2rad*lat1;
      		  double dlat = lat1r - lat0r;
      		  double dlon = lon1r - lon0r;
		  double y = sin(dlon) * cos(lat1r);
		  double x = cos(lat0r) * sin(lat1r) -
			sin(lat0r) * cos(lat1r) * cos(dlon);
		  double theta = atan2(y,x);
		  double baz = theta * 180 / M_PI;
		  return baz<0? 360 + baz: baz;
	    }
	    bool intersectLat(double inLat, double & outLon)
	    {
		  double b;
		  double lon0r = deg2rad*lon0;
      		  double lat0r = deg2rad*lat0;
      		  double lat3r = deg2rad*inLat;
		  double baz15 = (lon1>lon0) ? initialAzimuth()*deg2rad:(2*M_PI-initialAzimuth()*deg2rad);
		  double x = sin(baz15) * cos(lat0r)/cos(lat3r);
		  if(x==1)
			b = M_PI / 2; // one triangle
		  else if(x<1)
			b = (inLat<lat0) ? asin(x):(M_PI-asin(x)); // two triangles
		  else 
			return false;

		  double t1 = cos((lat0r-lat3r)/2)*cos((baz15+b)/2);
		  double t2 = sin((lat0r+lat3r)/2)*sin((baz15+b)/2);
		  double t3 = (atan2(t1,t2)>0) ? (atan2(t1,t2)*rad2deg) : (360+atan2(t1,t2)*rad2deg);
		  outLon = (lon0 > lon1) ? (lon0-2*t3) : (lon0+2*t3);
		  return (outLon-lon0)*(outLon-lon1)<0?true:false;
	    }
	    bool intersectLon(double inLon, double & outLat)
	    {
		  if((lon0-inLon)*(lon1-inLon)>0)
			return false;
		  double lon0r = deg2rad*lon0;
      		  double lat0r = deg2rad*lat0;
      		  double lon3r = deg2rad*inLon;
      		  double dlon = (lon3r>lon0r)?(lon3r-lon0r):(lon0r-lon3r);
		  double baz15 = (lon1>lon0) ? initialAzimuth()*deg2rad:(2*M_PI-initialAzimuth()*deg2rad);

		  double  x = -cos(baz15)*cos(dlon) + sin(baz15)*sin(dlon)*sin(lat0r);
		  double b = acos(x);
		  double t1 = cos(lat0r)*sin(baz15)*sin(dlon);
		  double t2 = cos(baz15)+cos(dlon)*cos(b);
		  double t3 = (atan2(t1,t2)>0) ? (atan2(t1,t2)*rad2deg) : (360+atan2(t1,t2)*rad2deg);
		  outLat = 90 - t3;
		  return true;
	    }
};

bool ge(double a1, double a2)
{
      return ((a1-a2) > -epsilon) ? true:false;
}
bool le(double a1, double a2)
{
      return ((a2-a1) > -epsilon) ? true:false;
}
void getData(std::string fileName, std::vector<std::vector<double>> & out)
{
      // open file
      std::ifstream in(fileName.c_str());
      if(!in)
	    std::cerr << "Cannot open the file: " << fileName << std::endl;
      std::string str;
      //load contend line by line
      while (std::getline(in,str))
      {
	    std::vector<double> row;
	    std::istringstream s(str);
	    double item;
	    while (s >> item)
		  row.push_back(item);
	    out.push_back(row);
      }
      in.close();
}

void mesh2D(std::vector<std::vector<double>> in, double dlon,double dlat, std::vector<double> & meshLon, std::vector<double> & meshLat)
{
      std::vector<double> lon,lat;
      for(int i=0;i<in.size();i++)
      {
	    lon.push_back(in[i][0]);
	    lon.push_back(in[i][2]);
	    lat.push_back(in[i][1]);
	    lat.push_back(in[i][3]);
      }
      std::sort(lon.begin(),lon.end());
      std::sort(lat.begin(),lat.end());
      int no = ceil ((lon.back()-lon.front()) / dlon )+ 1;
      for(int i=0;i<no+1;i++)
	    meshLon.push_back(lon.front()-0.5*dlon+i*dlon);

      no = ceil((lat.back()-lat.front())/dlat) + 1;
      for(int i=0;i<no+1;i++)
	    meshLat.push_back(lat.front()-0.5*dlat+i*dlat);
}	    
	    
void rayTrace(Rays ray, std::vector<double> meshLon, std::vector<double> meshLat, std::vector<std::vector<double>> & segRay)
{
      double xmesh1,xmesh2,ymesh1,ymesh2;
      double interlat1,interlat2;
      double interlon1,interlon2;
      double t;
      int lonsize = meshLon.size()-1;
      int latsize = meshLat.size()-1;
      for(int i=0;i<lonsize;i++)
      {
	    std::vector<double> row(latsize,0);
	    Rays tempRay;
	    xmesh1 = meshLon[i];
	    xmesh2 = meshLon[i+1];
	    interlat1 = ray.intersectLon(xmesh1,t)?t:1000;
	    interlat2 = ray.intersectLon(xmesh2,t)?t:1000;
	    for(int j=0;j<latsize;j++)
	    {
		  ymesh1 = meshLat[j];
		  ymesh2 = meshLat[j+1];
		  std::vector<double> intersections;
		  interlon1 = ray.intersectLat(ymesh1,t)?t:1000;
  		  interlon2 = ray.intersectLat(ymesh2,t)?t:1000;
		  if(ge(ray.lon0,xmesh1) && ge(xmesh2,ray.lon0)
			      && ge(ray.lat0,ymesh1) && ge(ymesh2,ray.lat0))
		  {
			intersections.push_back(ray.lon0);
			intersections.push_back(ray.lat0);

		  }
		  if(ge(ray.lon1,xmesh1) && ge(xmesh2,ray.lon1)
			      && ge(ray.lat1,ymesh1) && ge(ymesh2,ray.lat1))
		  {
			intersections.push_back(ray.lon1);
			intersections.push_back(ray.lat1);

		  }
		  if(ge(interlat1,ymesh1) && ge(ymesh2,interlat1))
		  {
		  	intersections.push_back(xmesh1);
	  		intersections.push_back(interlat1);
	    	  }
		  if(ge(interlat2,ymesh1) && ge(ymesh2,interlat2))
		  {
		  	intersections.push_back(xmesh2);
	  		intersections.push_back(interlat2);
	    	  }
		  if(ge(interlon1,xmesh1) && ge(xmesh2,interlon1))
	  	  {
		  	intersections.push_back(interlon1);
	  		intersections.push_back(ymesh1);
		  }
		  if(ge(interlon2,xmesh1) && ge(xmesh2,interlon2))
	  	  {
		  	intersections.push_back(interlon2);
	  		intersections.push_back(ymesh2);
		  }
		  if(intersections.size() == 4)
		  {
			tempRay.lon0 = intersections[0];
			tempRay.lat0 = intersections[1];
			tempRay.lon1 = intersections[2];
			tempRay.lat1 = intersections[3];
			row[j] = tempRay.length();

		  }

	    }
	    segRay.push_back(row);
      }
}


int main(int argc, char *argv[])
{	
      if (argc < 5)
      {
            std::cerr << "Usage: " << argv[0] << " <input_data_file> <output_results_file> <longitude_increment> <latitude_increment> [optional: predefined_mesh_file]\n";
            return 1; // Return error code
      }
      std::vector<std::vector<double>> data;
      std::vector<double> meshLon, meshLat;
      Rays ray;
      getData(argv[1],data);
      std::ofstream outFile(argv[2]);
      double dlon = atof(argv[3]);
      double dlat = atof(argv[4]);
      if(argc == 6)
      {
	    std::vector<std::vector<double>> mtemp;
	    getData(argv[5],mtemp);
	    for (int i=0;i<mtemp.size();i++)
	    {
		  meshLon.push_back(mtemp[i][0]);
		  meshLat.push_back(mtemp[i][1]);
	    }
	    std::sort(meshLon.begin(),meshLon.end());
	    meshLon.erase(unique(meshLon.begin(),meshLon.end()),meshLon.end() );
	    std::sort(meshLat.begin(),meshLat.end());
	    meshLat.erase(unique(meshLat.begin(),meshLat.end()),meshLat.end() );
      }
      else if(argc == 5) 
      {
	    mesh2D(data,dlon,dlat,meshLon,meshLat);
      }
      for(int i=0;i<meshLon.size()-1;++i)
      {
	    for(int j=0;j<meshLat.size()-1;++j)
		  std::cout << meshLon[i] << " " << meshLat[j] << std::endl;
      }
      for(int i=0;i<data.size();i++)
      {
	    std::vector<std::vector<double>> segRay;
	    ray.lon0 = data[i][0];
	    ray.lat0 = data[i][1];
	    ray.lon1 = data[i][2];
	    ray.lat1 = data[i][3];
	    rayTrace(ray,meshLon,meshLat,segRay);
	    double l = 0.0;
	    double ar = 0;
	    for(int r=0;r<segRay.size();r++)
	    {
		  for(int c=0;c<segRay[0].size();c++)
			l += segRay[r][c];
	    }
		  for (auto item : data[i])
			outFile << item << " ";

		  for(int r=0;r<segRay.size();r++)
		  {
			for(int c=0;c<segRay[0].size();c++)
			      outFile << segRay[r][c] << " ";
		  }
		  outFile << "\n";
      }
}

