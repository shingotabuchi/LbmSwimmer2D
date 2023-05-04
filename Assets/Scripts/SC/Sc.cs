using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Sc : MonoBehaviour
{
    public int nx=128;
    public int ny=128;
    public Image plotImage;
    Texture2D plotTexture;
    RenderTexture renderTexture;
    public float tau=1.0f;
    public float rhol=1.95f;
    public float rhog=0.15f;
    public int radius=20;
    public float g=-5.0f;
    //Arrays
    int npop=9;
    ComputeBuffer rho,u1,u2,f;
    // Color[] pixels;
    public ComputeShader compute;
    int init,plotDensity,collisionStreaming;
    private void OnValidate() {
        compute.SetFloat("tau",tau);
        compute.SetFloat("rhol",rhol);
        compute.SetFloat("rhog",rhog);
        compute.SetFloat("g",g);
    }
    private void OnDestroy() {
        rho.Dispose();
        u1.Dispose();
        u2.Dispose();
        f.Dispose();
    }    
    void Start()
    {
        plotTexture = new Texture2D(nx,ny);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,nx,ny),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(1080,1080);
        renderTexture = new RenderTexture(nx,ny,24);
        renderTexture.enableRandomWrite = true;

        rho = new ComputeBuffer(nx*ny,sizeof(float));
        u1 = new ComputeBuffer(nx*ny,sizeof(float));
        u2 = new ComputeBuffer(nx*ny,sizeof(float));
        f = new ComputeBuffer(nx*ny*npop*2,sizeof(float));
        compute.SetInt("npop",npop);
        compute.SetInt("nx",nx);
        compute.SetInt("ny",ny);
        compute.SetInt("nxy9",nx*ny*9);
        compute.SetInt("radius",radius);
        OnValidate();

        init = compute.FindKernel("init");
        compute.SetBuffer(init,"rho",rho);
        compute.SetBuffer(init,"u1",u1);
        compute.SetBuffer(init,"u2",u2);
        compute.SetBuffer(init,"f",f);

        plotDensity = compute.FindKernel("plotDensity");
        compute.SetBuffer(plotDensity,"f",f);
        compute.SetBuffer(plotDensity,"rho",rho);
        compute.SetTexture(plotDensity,"renderTexture",renderTexture);

        collisionStreaming = compute.FindKernel("collisionStreaming");
        compute.SetBuffer(collisionStreaming,"rho",rho);
        compute.SetBuffer(collisionStreaming,"u1",u1);
        compute.SetBuffer(collisionStreaming,"u2",u2);
        compute.SetBuffer(collisionStreaming,"f",f);

        compute.Dispatch(init,(nx*ny+63)/64,1,1);
    }

    private void FixedUpdate() 
    {
        compute.Dispatch(plotDensity,(nx*ny+63)/64,1,1);
        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();

        compute.Dispatch(collisionStreaming,(nx+7)/8,(ny+7)/8,1);
    }
}
