using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class AbbasCs : MonoBehaviour
{
    public int DIM;
    public Image plotImage;
    Texture2D plotTexture;

	public float minPhase,maxPhase;
	public float minRho,maxRho;
	public float gravity;

	int tf   = 10000;
	int step;
	
	int X0;
	int Y0;

	int[] ex = new int[9]{0, 1, 0,-1, 0, 1,-1,-1, 1};
	int[] ey = new int[9]{0, 0, 1, 0,-1, 1, 1,-1,-1};
	float[] Wa = new float[9]{16f/36f,4f/36f, 4f/36f, 4f/36f, 4f/36f, 1f/36f, 1f/36f, 1f/36f, 1f/36f};
	float R;
	float Rhol = 0.001f;
	float Rhoh = 1f;
	float dRho3;
	float Sigma = 0.01f;
	float tau =  0.3f + 0.5f;
	float s8;
	float W = 4;
	float Beta,K;
	float M = 0.02f;
	float w_c;
	float[] h,g,C,P,mu,DcDx,DcDy,Rho,Ux,Uy,ni,nj;
	float[] Gamma = new float[9];
	float[] Ga_Wa = new float[9];
	float[] heq = new float[9];
	float[] geq = new float[9];
	float[] hlp = new float[9];
	float[] eF = new float[9];
	float tmp,Ri,Fx,Fy;
	int DIMSqrd9;
    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM,DIM);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM,DIM),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(1080,1080);

		X0 = DIM/2 + 1;
		Y0 = DIM/2 + 1;
		R = (float)(DIM)/8f;
		step = tf/10;
		dRho3 = (Rhoh - Rhol)/3f;
		s8 = 1f/tau;
		Beta = 12f * Sigma/W;
		K    = 1.5f * Sigma*W;
		w_c = 1f/(0.5f + 3f*M);
		DIMSqrd9 = DIM*DIM*9;
		h = new float[DIM*DIM*9*2];
		g = new float[DIM*DIM*9*2];
		C = new float[DIM*DIM];
		P = new float[DIM*DIM];
		mu = new float[DIM*DIM];
		DcDx = new float[DIM*DIM];
		DcDy = new float[DIM*DIM];
		Rho = new float[DIM*DIM];
		Ux = new float[DIM*DIM];
		Uy = new float[DIM*DIM];
		ni = new float[DIM*DIM];
		nj = new float[DIM*DIM];
		Initialize_distributions();
		// PlotPhase();
		PlotRho();
	}

	void Initialize_distributions()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				int index = ArrayIndex(i,j);
				P[index] = 0f;
				Ux[index] = 0f;
				Uy[index] = 0f;
				// if(Mathf.Abs(i-(X0-0.5f)) < R && Mathf.Abs(j-(Y0-0.5f)) < R)
				// {
					// C[index] = Random.Range(-0.01f,0.01f);
				// }
				// else C[index] = 1;
				Ri = Mathf.Sqrt( (i-(X0-0.5f))*(i-(X0-0.5f)) + (j-(Y0-0.5f))*(j-(Y0-0.5f)) );
				C[index] = 0.5f - 0.5f * (float)System.Math.Tanh(2f*(R-Ri)/W);
			}
		}

		// ChemPotAndGradients();
		// normal_FD();

		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				int index = ArrayIndex(i,j);

				DcDx[index] = (C[ArrayIndex(i+1,j  )] - C[ArrayIndex(i-1,j  )])/3f + ( C[ArrayIndex(i+1,j-1)] + C[ArrayIndex(i+1,j+1)] - C[ArrayIndex(i-1,j-1)] - C[ArrayIndex(i-1,j+1)])/12f;
				DcDy[index] = (C[ArrayIndex(i  ,j+1)] - C[ArrayIndex(i  ,j-1)])/3f + ( C[ArrayIndex(i-1,j+1)] + C[ArrayIndex(i+1,j+1)] - C[ArrayIndex(i-1,j-1)] - C[ArrayIndex(i+1,j-1)])/12f;
				float D2C = ( C[ArrayIndex(i-1,j-1)]+C[ArrayIndex(i+1,j-1)]+C[ArrayIndex(i-1,j+1)]+C[ArrayIndex(i+1,j+1)]+4*(C[ArrayIndex(i  ,j-1)]+C[ArrayIndex(i-1,j  )]+C[ArrayIndex(i+1,j  )]+C[ArrayIndex(i  ,j+1)]) - 20*C[ArrayIndex(i,j)] )/6f;
				mu[index] = 4f*Beta * C[index] * (C[index]-1f) * (C[index]-0.5f) - K*D2C;

				tmp = Mathf.Sqrt(DcDx[index]*DcDx[index] + DcDy[index]*DcDy[index] + 1e-32f);
				ni[index] = DcDx[index] / tmp;
				nj[index] = DcDy[index] / tmp;

				Rho[index] = Rhol + C[index] * (Rhoh - Rhol);
				P[index] = P[index] - C[index] * Sigma/R /(Rho[index]/3);
				// Equilibrium_new( Ux[index], Uy[index] );
				float u = Ux[index];
				float v = Uy[index];
				float u2 = u*u + v*v;
				for (int k = 0; k < 9; k++)
				{
					float eU = ex[k] * u  + ey[k] * v;
					Ga_Wa[k] = Wa[k] * ( eU*(3f + 4.5f*eU) - 1.5f*u2 );
					Gamma[k] = Ga_Wa[k] + Wa[k];
					eF[k]  = ( 1f - 4f*(C[index]-0.5f)*(C[index]-0.5f) )/W * ( ex[k]*ni[index] + ey[k]*nj[index] );
					hlp[k] = Wa[k] * eF[k];
					h[index*9 + k] = C[index]*Gamma[k] - 0.5f * hlp[k];
					g[index*9 + k] = P[index] * Wa[k] + Ga_Wa[k];
				}
			}
		}
	}

	void PlotPhase()
	{
		Color[] pixels = new Color[DIM*DIM];
		for (int i = 0; i < DIM*DIM; i++)
		{
			pixels[i] = ColorMap(C[i]-minPhase,maxPhase - minPhase);
		}

		plotTexture.SetPixels(pixels);
		plotTexture.Apply();
	}

	void PlotRho()
	{
		Color[] pixels = new Color[DIM*DIM];
		for (int i = 0; i < DIM*DIM; i++)
		{
			pixels[i] = ColorMap(Rho[i]-minRho,maxRho - minRho);
		}

		plotTexture.SetPixels(pixels);
		plotTexture.Apply();
	}

	// void Equilibrium_new(float u, float v)
	// {
	// 	float u2 = u*u + v*v;

	// 	for (int i = 0; i < 9; i++)
	// 	{
	// 		float eU = ex[i] * u  + ey[i] * v;
	// 		Ga_Wa[i] = Wa[i] * ( eU*(3f + 4.5f*eU) - 1.5f*u2 );
	// 	}
	// }

	// void ChemPotAndGradients()
	// {
	// 	float D2C;
	// 	for (int i = 0; i < DIM; i++)
	// 	{
	// 		for (int j = 0; j < DIM; j++)
	// 		{
	// 			int index = ArrayIndex(i,j);
	// 			DcDx[index] = (C[ArrayIndex(i+1,j  )] - C[ArrayIndex(i-1,j  )])/3f + ( C[ArrayIndex(i+1,j-1)] + C[ArrayIndex(i+1,j+1)] - C[ArrayIndex(i-1,j-1)] - C[ArrayIndex(i-1,j+1)])/12f;
	// 			DcDy[index] = (C[ArrayIndex(i  ,j+1)] - C[ArrayIndex(i  ,j-1)])/3f + ( C[ArrayIndex(i-1,j+1)] + C[ArrayIndex(i+1,j+1)] - C[ArrayIndex(i-1,j-1)] - C[ArrayIndex(i+1,j-1)])/12f;
	// 			D2C = ( C[ArrayIndex(i-1,j-1)]+C[ArrayIndex(i+1,j-1)]+C[ArrayIndex(i-1,j+1)]+C[ArrayIndex(i+1,j+1)]+4*(C[ArrayIndex(i  ,j-1)]+C[ArrayIndex(i-1,j  )]+C[ArrayIndex(i+1,j  )]+C[ArrayIndex(i  ,j+1)]) - 20*C[ArrayIndex(i,j)] )/6f;
	// 			mu[index] = 4f*Beta * C[index] * (C[index]-1f) * (C[index]-0.5f) - K*D2C;
	// 		}
	// 	}
	// }

	int ArrayIndex(int x,int y)
	{
		return (x+DIM)%DIM + ((y+DIM)%DIM)*DIM;
	}

    // Update is called once per frame
    void FixedUpdate()
    {
        Collision_h_g();
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				int index = ArrayIndex(i,j);
				for (int k = 0; k < 9; k++)
				{
					h[DIMSqrd9 + index * 9 + k] = h[ArrayIndex(i-ex[k],j-ey[k])*9 + k];
					g[DIMSqrd9 + index * 9 + k] = g[ArrayIndex(i-ex[k],j-ey[k])*9 + k];
				}
			}
		}
		Macroscopic_Properties_h();
		Macroscopic_Properties_g();
		PlotRho();
    }

	void Macroscopic_Properties_h()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				int index = ArrayIndex(i,j);
				C[index] = 0;
				for (int k = 0; k < 9; k++)
				{
					h[index*9 + k] = h[DIMSqrd9 + index * 9 + k];
					C[index] += h[index*9 + k];
				}
				Rho[index] = Rhol + C[index] * (Rhoh - Rhol);
			}
		}
	}

	void Macroscopic_Properties_g()
	{
		// Boundary_Conditions_C();
		// ChemPotAndGradients();
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				int index = ArrayIndex(i,j);
				DcDx[index] = (C[ArrayIndex(i+1,j  )] - C[ArrayIndex(i-1,j  )])/3f + ( C[ArrayIndex(i+1,j-1)] + C[ArrayIndex(i+1,j+1)] - C[ArrayIndex(i-1,j-1)] - C[ArrayIndex(i-1,j+1)])/12f;
				DcDy[index] = (C[ArrayIndex(i  ,j+1)] - C[ArrayIndex(i  ,j-1)])/3f + ( C[ArrayIndex(i-1,j+1)] + C[ArrayIndex(i+1,j+1)] - C[ArrayIndex(i-1,j-1)] - C[ArrayIndex(i+1,j-1)])/12f;
				float D2C = ( C[ArrayIndex(i-1,j-1)]+C[ArrayIndex(i+1,j-1)]+C[ArrayIndex(i-1,j+1)]+C[ArrayIndex(i+1,j+1)]+4*(C[ArrayIndex(i  ,j-1)]+C[ArrayIndex(i-1,j  )]+C[ArrayIndex(i+1,j  )]+C[ArrayIndex(i  ,j+1)]) - 20*C[ArrayIndex(i,j)] )/6f;
				mu[index] = 4f*Beta * C[index] * (C[index]-1f) * (C[index]-0.5f) - K*D2C;

				P[index] = 0;
				for (int k = 0; k < 9; k++)
				{
					g[index*9 + k] = g[DIMSqrd9 + index * 9 + k];
					P[index] += g[index*9 + k];
				}

				float FpX = - P[index] * dRho3 * DcDx[index];
				float FpY = - P[index] * dRho3 * DcDy[index];

				float u = Ux[index];
				float v = Uy[index];
				// Equilibrium_new( Ux[index], Uy[index] );
				float u2 = u*u + v*v;
				float sxx = 0f;
				float sxy = 0f;
				float syy = 0f;
				for (int k = 0; k < 9; k++)
				{
					float eU = ex[k] * u  + ey[k] * v;
					Ga_Wa[k] = Wa[k] * ( eU*(3f + 4.5f*eU) - 1.5f*u2 );
					geq[k] = P[index] * Wa[k] + Ga_Wa[k];

					sxx +=  (g[index*9 + k]-geq[k])*ex[k]*ex[k];
					sxy +=  (g[index*9 + k]-geq[k])*ex[k]*ey[k];
					sxy +=  (g[index*9 + k]-geq[k])*ey[k]*ey[k];
				}
				// float FmX = 0;
				// float FmY = 0;
				float FmX = (0.5f-tau)/tau * (sxx*DcDx[index]+sxy*DcDy[index]) * (Rhoh-Rhol);
				float FmY = (0.5f-tau)/tau * (sxy*DcDx[index]+syy*DcDy[index]) * (Rhoh-Rhol);
				// Calculate_Viscous_Force( index, DcDx[index], DcDy[index],ref FmX,ref FmY );

				Fx = mu[index] * DcDx[index] + FpX + FmX;
				Fy = mu[index] * DcDy[index] + FpY + FmY - gravity * Rho[index];

				Ux[index] = g[index*9 + 1]-g[index*9 + 3]+g[index*9 + 5]-g[index*9 + 6]-g[index*9 + 7]+g[index*9 + 8] + 0.5f*Fx/Rho[index];
				Uy[index] = g[index*9 + 2]-g[index*9 + 4]+g[index*9 + 5]+g[index*9 + 6]-g[index*9 + 7]-g[index*9 + 8] + 0.5f*Fy/Rho[index];
			}
		}
	}

	// void Boundary_Conditions_C()
	// {
	// 	for (int i = 0; i < DIM; i++)
	// 	{
	// 		C[i + 0*DIM] = C[i + (DIM-2)*DIM];
	// 		C[i + (DIM-1)*DIM] = C[i + 1*DIM];
	// 		C[0 + i*DIM] = C[DIM-2 + i*DIM];
	// 		C[DIM-1 + i*DIM] = C[1 + i*DIM];
	// 	}
	// }

	float[] Streaming(float[] f)
	{
		float[] fnew = new float[9*DIM*DIM];
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				int index = ArrayIndex(i,j);
				for (int k = 0; k < 9; k++)
				{
					fnew[index * 9 + k] = f[ArrayIndex(i-ex[k],j-ey[k])*9 + k];
				}
			}
		}
		return fnew;
	}

	// void Boundary_Conditions_f(float[] f)
	// {
	// 	for (int i = 0; i < DIM; i++)
	// 	{
	// 		for (int k = 0; k < 9; k++)
	// 		{
	// 			f[k + (0 + i*DIM)*9] = f[k + (DIM-2 + i*DIM)*9];
	// 			f[k + (DIM-1 + i*DIM)*9] = f[k + (1 + i*DIM)*9];

	// 			f[k + (i + 0*DIM)*9] = f[k + (i + (DIM-2)*DIM)*9];
	// 			f[k + (i + (DIM-1)*DIM)*9] = f[k + (i + 1*DIM)*9];
	// 		}
	// 	}
	// }

	void Collision_h_g()
	{
		float FpX, FpY, FmX, FmY;
		FpX = 0;
		FpY = 0;
		FmX = 0;
		FmY = 0;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				int index = ArrayIndex(i,j);
				tmp = Mathf.Sqrt(DcDx[index]*DcDx[index] + DcDy[index]*DcDy[index] + 1e-32f);
				ni[index] = DcDx[index] / tmp;
				nj[index] = DcDy[index] / tmp;

				// Equilibrium_new( Ux[index], Uy[index] );
				float u = Ux[index];
				float v = Uy[index];
				float u2 = u*u + v*v;
				float sxx = 0f;
				float sxy = 0f;
				float syy = 0f;
				for (int k = 0; k < 9; k++)
				{
					float eU = ex[k] * u  + ey[k] * v;
					Ga_Wa[k] = Wa[k] * ( eU*(3f + 4.5f*eU) - 1.5f*u2 );
					Gamma[k] = Ga_Wa[k] + Wa[k];

					eF[k]  = ( 1f - 4f*(C[index]-0.5f)*(C[index]-0.5f) )/W * ( ex[k]*ni[index] + ey[k]*nj[index] );

					hlp[k] = Wa[k] * eF[k];

					heq[k] = C[index]*Gamma[k] - 0.5f * hlp[k];

					h[index*9 +k] = h[index*9 +k] * (1f-w_c) + heq[k] * w_c + hlp[k];

					FpX = - P[index] * dRho3 * DcDx[index];
					FpY = - P[index] * dRho3 * DcDy[index];

					geq[k] = P[index] * Wa[k] + Ga_Wa[k];
					sxx +=  (g[index*9 + k]-geq[k])*ex[k]*ex[k];
					sxy +=  (g[index*9 + k]-geq[k])*ex[k]*ey[k];
					sxy +=  (g[index*9 + k]-geq[k])*ey[k]*ey[k];
					
				}
				// for (int k = 0; k < 9; k++)
				// {
				// 	geq[k] = P[index] * Wa[k] + Ga_Wa[k];
				// }
				FmX = (0.5f-tau)/tau * (sxx*DcDx[index]+sxy*DcDy[index]) * (Rhoh-Rhol);
				FmY = (0.5f-tau)/tau * (sxy*DcDx[index]+syy*DcDy[index]) * (Rhoh-Rhol);
				Fx = mu[index] * DcDx[index] + FpX + FmX;
				Fy = mu[index] * DcDy[index] + FpY + FmY - gravity * Rho[index];
				for (int k = 0; k < 9; k++)
				{
					// Calculate_Viscous_Force(index, DcDx[index], DcDy[index], ref FmX,ref FmY );
					eF[k] = ex[k] * Fx + ey[k] * Fy;

					hlp[k] = 3f * Wa[k] * eF[k] / Rho[index];

					geq[k] = P[index] * Wa[k] + Ga_Wa[k] - 0.5f * hlp[k];

					g[index*9 +k] = g[index*9 +k] * (1f-s8) + geq[k] * s8 + hlp[k];
				}
			}
		}
	}

	// void Calculate_Viscous_Force(int index, float _DcDx,float _DcDy,ref float FmX,ref float FmY )
	// {
	// 	Calculate_Viscous_Force_BGK(index,_DcDx, _DcDy,ref FmX,ref FmY );
	// }

	// void Calculate_Viscous_Force_BGK(int index, float _DcDx,float _DcDy,ref float FmX,ref float FmY )
	// {
	// 	float sxx = 0f, sxy = 0f, syy = 0f;
	// 	Calculate_Stress_Tensor_BGK(index, ref sxx,ref sxy,ref syy );
	// 	FmX = (0.5f-tau)/tau * (sxx*_DcDx+sxy*_DcDy) * (Rhoh-Rhol);
	// 	FmY = (0.5f-tau)/tau * (sxy*_DcDx+syy*_DcDy) * (Rhoh-Rhol);
	// }

	// void Calculate_Stress_Tensor_BGK(int index, ref float sxx,ref float sxy,ref float syy )
	// {
	// 	sxx = 0f;
	// 	sxy = 0f;
	// 	syy = 0f;
	// 	for (int i = 1; i < 9; i++)
	// 	{
	// 		sxx +=  (g[index*9 + i]-geq[i])*ex[i]*ex[i];
	// 		sxy +=  (g[index*9 + i]-geq[i])*ex[i]*ey[i];
	// 		sxy +=  (g[index*9 + i]-geq[i])*ey[i]*ey[i];
	// 	}
	// }

	Color[] colorsOfMap = 
	new Color[]{
		new Color(0,0,0,1),
		new Color(0,0,1,1),
		new Color(0,1,1,1),
		new Color(0,1,0,1),
		new Color(1,1,0,1),
		new Color(1,0,0,1),
		new Color(1,1,1,1)
	};

	float colorPerc = 1.0f / (7.0f-1.0f);

	Color ColorMap(float val, float maxVal)
	{
		if(val > maxVal) val = maxVal;
		float valPerc = val / maxVal;// value%
		int blockIdx = (int)(valPerc / colorPerc);// Idx of 
		float valPercResidual = valPerc - (blockIdx*colorPerc);//remove the part represented of block 
		float percOfColor = valPercResidual / colorPerc;// % of color of this block that will be filled
		Color cTarget = colorsOfMap[blockIdx];
		float deltaR = 0;
		float deltaG = 0;
		float deltaB = 0;

		if(blockIdx != 6)
		{
			Color cNext = colorsOfMap[blockIdx + 1];
			deltaR =cNext.r - cTarget.r;
			deltaG =cNext.g - cTarget.g;
			deltaB =cNext.b - cTarget.b;
		}

		float R = cTarget.r + (deltaR * percOfColor);
		float G = cTarget.g + (deltaG * percOfColor);
		float B = cTarget.b + (deltaB * percOfColor);

		return new Color(R,G,B,1.0f);
	}
}
