#include "Altitude.h"
#include<math.h>
#include<iostream>

Altitude::Altitude(double h) : m_height(static_cast<float>(h))
{
	m_rho	= m_rhoSeaLevel * exp(-0.0296f * m_height / 304.8f);
	m_temp	= (15.0f - 0.001981f * m_height / 304.8f) + 273.16f;
	m_a		= pow(1.4f * R * m_temp, 0.5f);
	m_p		= m_pSeaLevel * pow((1 - 6.876e-6f * m_height * 3.2808f), 5.265f);
	m_vis = m_visSeaLevel * ((m_tSeaLevel + 110.0f) / (m_temp + 110.0f)) * pow(m_temp / m_tSeaLevel, 1.5f);
}


void Altitude::setAltitude(double h)
{
	m_height	= static_cast<float>(h);
	m_rho		= m_rhoSeaLevel * exp(-0.0296f * m_height / 304.8f);
	m_temp		= (15.0f - 0.001981f * m_height / 304.8f) + 273.16f;
	m_a			= pow(1.4f * R * m_temp, 0.5f);
	m_p			= m_pSeaLevel * pow((1 - 6.876e-6f * m_height * 3.2808f), 5.265f);
	m_vis		= m_visSeaLevel * ((m_tSeaLevel + 110.0f) / (m_temp + 110.0f)) * pow(m_temp / m_tSeaLevel, 1.5f);
}

float Altitude::getrho() const
{
	return m_rho;
}

float Altitude::getSOS() const
{
	return m_a;
}

float Altitude::getTemp() const
{
	return m_temp;
}

float Altitude::getVis() const
{
	return m_vis;
}

float Altitude::getPres() const
{
	return m_p;
}

