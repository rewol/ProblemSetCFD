#ifndef ALTITUDE
#define ALTITUDE

class Altitude
{
friend class Field2D;

public:
	Altitude(double);
	void setAltitude(double);
	float getrho() const;
	float getSOS() const;
	float getTemp() const;
	float getVis() const;
	float getPres() const;

private:
	float m_height;
	float m_rho;
	float m_temp;
	float m_a;
	float m_p;
	float m_vis;
	const float m_rhoSeaLevel		= 1.225f;		// In kg/m3
	const float m_tSeaLevel			= 288.16f;		// In Kelvin
	const float R					= 287.0f;		// J/(kg.K)
	const float m_visSeaLevel		= 1.7894e-5f;	// kg/(m.s)
	const float m_pSeaLevel			= 101301.0f;	// N/m2

};

#endif