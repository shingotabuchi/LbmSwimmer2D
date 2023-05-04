using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Sc : MonoBehaviour
{
    public int nx=128;
    public int ny=128;
    public float maxRho = 1f;
    public Image plotImage;
    Texture2D plotTexture;
    int[] cx=new int[]{0,1,0,-1,0,1,-1,-1,1};
    int[] cy=new int[]{0,0,1,0,-1,1,1,-1,-1};
    float[] weights=new float[]{4.0f/9.0f,1.0f/9.0f,1.0f/9.0f,1.0f/9.0f,1.0f/9.0f,1.0f/36.0f,1.0f/36.0f,1.0f/36.0f,1.0f/36.0f};
    float tau=1.0f;
    float rhol=1.95f;
    float rhog=0.15f;
    int radius=20;
    float g=-5.0f;
    //Arrays
    int npop=9;
    float[] rho,u1,u2,f,f2,feq;
    Color[] pixels;
    public int stepCount = 0;
    void Start()
    {
        plotTexture = new Texture2D(nx,ny);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,nx,ny),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(1080,1080);
	
        rho = new float[nx*ny];
        u1 = new float[nx*ny];
        u2 = new float[nx*ny];
        f = new float[nx*ny*npop];
        f2 = new float[nx*ny*npop];
        feq = new float[npop];
        pixels = new Color[nx*ny];

        //Initialization
        for (int i=0; i < nx*ny; i++)
        {	
            // if ((i/nx - ny/2.0f)*(i/nx - ny/2.0f) + (i%nx - nx/2.0f)*(i%nx - nx/2.0f) <= radius * radius)
            if ((i/nx - ny/2.0f)*(i/nx - ny/2.0f)<= radius * radius&&(i%nx - nx/2.0f)*(i%nx - nx/2.0f) <= radius * radius)
            {
                rho[i]=rhol;
            }
            else 
                rho[i]=rhog;

            float dense,v1,v2;
            
            dense=rho[i];
            v1=v2=u1[i]=u2[i]=0.0f;
            float usq = v1*v1 + v2*v2;
            feq[0] = 4.0f/9.0f * dense * (1.0f - 1.5f * usq); 
            feq[1] = 1.0f/9.0f * dense * (1.0f + 3*v1 + 4.5f*v1*v1 - 1.5f*usq); 
            feq[2] = 1.0f/9.0f * dense * (1.0f + 3*v2 + 4.5f*v2*v2 - 1.5f*usq); 
            feq[3] = 1.0f/9.0f * dense * (1.0f - 3*v1 + 4.5f*v1*v1 - 1.5f*usq); 
            feq[4] = 1.0f/9.0f * dense * (1.0f - 3*v2 + 4.5f*v2*v2 - 1.5f*usq); 
            feq[5] = 1.0f/36.0f * dense * (1.0f + 3*(v1 + v2) + 4.5f*(v1 + v2)*(v1 + v2) - 1.5f*usq); 
            feq[6] = 1.0f/36.0f * dense * (1.0f + 3*(-v1 + v2) + 4.5f*(-v1 + v2)*(-v1 + v2) - 1.5f*usq);
            feq[7] = 1.0f/36.0f * dense * (1.0f + 3*(-v1 - v2) + 4.5f*(v1 + v2)*(v1 + v2) - 1.5f*usq); 
            feq[8] = 1.0f/36.0f * dense * (1.0f + 3*(v1 - v2) + 4.5f*(v1 - v2)*(v1 -v2) - 1.5f*usq); 
            for (int k=0; k<npop; k++) {
                f[9*i+k]=feq[k];
                f2[9*i+k]=feq[k];
            }
        }
    }

    private void FixedUpdate() 
    {
        if(stepCount == 30000)
        {
            print("dunzo");
            return;
        }
        stepCount++;
        //Calculation of the density field
        
		for (int i=0; i<nx*ny; i++) 
		{
			rho[i]=0; 
			for (int k=0; k<9; k++ )
			{			
				rho[i]+=f[9*i+k]; 
			}		
            pixels[i] = ColorMap(rho[i],maxRho);
		}
        plotTexture.SetPixels(pixels);
        plotTexture.Apply();
        //Collision and streaming
		for (int iY=0; iY<ny; iY++) 
        {
		    for(int iX=0;iX<nx;iX++)
		    {
		    
				int i=iY*nx+iX;
				float dense,v1,v2;
		       	
				dense=rho[i];

				float fx=0.0f;
				float fy=0.0f;

				for(int k=0;k<9;k++)
				{
					int iX2=(iX+cx[k]+nx) % nx; 
					int iY2=(iY+cy[k]+ny) % ny;
					fx+=weights[k]*cx[k]*(1.0f-Mathf.Exp(-rho[nx*iY2+iX2]));
					fy+=weights[k]*cy[k]*(1.0f-Mathf.Exp(-rho[nx*iY2+iX2]));
				}
			
				fx=-g*(1.0f-Mathf.Exp(-rho[i]))*fx;
				fy=-g*(1.0f-Mathf.Exp(-rho[i]))*fy;
			
				v1=u1[i]=(f[9*i+1]-f[9*i+3]+f[9*i+5]-f[9*i+6]-f[9*i+7]+f[9*i+8])/dense+fx/(2.0f*dense); 
				v2=u2[i]=(f[9*i+2]-f[9*i+4]+f[9*i+5]+f[9*i+6]-f[9*i+7]-f[9*i+8])/dense+fy/(2.0f*dense); 
			
				float[] fpop = new float[9];
				for(int k=0;k<9;k++)
					fpop[k]=weights[k]*(1-0.5f/tau)*((3*(cx[k]-v1)+9*cx[k]*(cx[k]*v1+cy[k]*v2))*fx
		               +(3*(cy[k]-v2)+9*cy[k]*(cx[k]*v1+cy[k]*v2))*fy);
			
				float usq = v1*v1 + v2*v2;	
			
				feq[0] = 4.0f/9.0f * dense * (1.0f - 1.5f * usq); 
				feq[1] = 1.0f/9.0f * dense * (1.0f + 3*v1 + 4.5f*v1*v1 - 1.5f*usq); 
				feq[2] = 1.0f/9.0f * dense * (1.0f + 3*v2 + 4.5f*v2*v2 - 1.5f*usq); 
				feq[3] = 1.0f/9.0f * dense * (1.0f - 3*v1 + 4.5f*v1*v1 - 1.5f*usq); 
				feq[4] = 1.0f/9.0f * dense * (1.0f - 3*v2 + 4.5f*v2*v2 - 1.5f*usq); 
				feq[5] = 1.0f/36.0f * dense * (1.0f + 3*(v1 + v2) + 4.5f*(v1 + v2)*(v1 + v2) - 1.5f*usq); 
				feq[6] = 1.0f/36.0f * dense * (1.0f + 3*(-v1 + v2) + 4.5f*(-v1 + v2)*(-v1 + v2) - 1.5f*usq);
				feq[7] = 1.0f/36.0f * dense * (1.0f + 3*(-v1 - v2) + 4.5f*(v1 + v2)*(v1 + v2) - 1.5f*usq);
				feq[8] = 1.0f/36.0f * dense * (1.0f + 3*(v1 - v2) + 4.5f*(v1 - v2)*(v1 -v2) - 1.5f*usq);
			
				for(int k=0; k<9; k++) 
				{  
					int iX2=(iX+cx[k]+nx) % nx; 
					int iY2=(iY+cy[k]+nx) % nx;
					f[9*i+k]+=-1.0f/tau*(f[9*i+k]-feq[k])+fpop[k]; 
					f2[9*(nx*iY2+iX2)+k]=f[9*i+k]; 
				}  
			}
        }

        f = (float[])f2.Clone();
    }

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
