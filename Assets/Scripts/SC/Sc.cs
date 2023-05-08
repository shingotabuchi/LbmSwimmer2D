using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.VFX;
using UnityEngine.UI;

public class Sc : MonoBehaviour
{
    public enum InitMode
    {
        Random,
        Circle,
        Liquid,
        Gas,
        Droplet,
    }
    public enum PlotMode
    {
        Density,
        Temp,
        Speed,
    }
    public int nx=128;
    public int ny=128;
    public Image plotImage;
    Texture2D plotTexture;
    RenderTexture renderTexture;
    // RenderTexture renderTextureOne;
    // public float tau=1.0f;
    public float vl;
    public float vg;
    public float taut=1.0f;
    public float rhol=1.95f;
    public float rhog=0.15f;
    public float maxrho=1.95f;
    public float minrho=0.15f;
    public float maxtemp=1.95f;
    public float mintemp=0.15f;
    public int radius=20;
    public float gs=-0.27f;
    public float grav=-0.00003f;
    public float initTr = 0.86f;
    public float cv = 6f;
    float aveRho;
    public int loopCount;
    public bool periodic;
    public InitMode initMode;
    public PlotMode plotMode;
    public float minspeed,maxspeed;
    public float walltempr = 0.9f;
    public float wallrho;
    public VisualEffect vfx;
    GraphicsBuffer velocityGraphicsBuffer;
    //Arrays
    int npop=9;
    ComputeBuffer rhot,uv,f,g;
    // Color[] pixels;
    public ComputeShader compute;
    int init,collisionStreaming,bouncebackBoundary,oneLoop,plotPass,collisionStreamingG,plotSpeed;
    private void OnValidate() {
        compute.SetFloat("taut",taut);
        compute.SetFloat("vg",vg);
        compute.SetFloat("vl",vl);
        compute.SetFloat("rhol",rhol);
        compute.SetFloat("rhog",rhog);
        compute.SetFloat("maxrho",maxrho);
        compute.SetFloat("minrho",minrho);
        compute.SetFloat("maxtemp",maxtemp);
        compute.SetFloat("mintemp",mintemp);
        compute.SetFloat("maxspeed",maxspeed);
        compute.SetFloat("maxspeed",maxspeed);
        compute.SetFloat("wallrho",wallrho);
        compute.SetFloat("walltempr",walltempr);
        compute.SetFloat("gs",gs);
        compute.SetFloat("grav",grav);
        compute.SetBool("periodic",periodic);
        compute.SetInt("plotMode",(int)plotMode);
        compute.SetFloat("cv",cv);
    }
    private void OnDestroy() {
        rhot.Dispose();
        uv.Dispose();
        f.Dispose();
        g.Dispose();
    }    
    void Start()
    {
        // plotTexture = new Texture2D(1,1);
        plotTexture = new Texture2D(nx,ny);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,nx,ny),UnityEngine.Vector2.zero);
        // plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,1,1),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(1080,1080);
        renderTexture = new RenderTexture(nx,ny,24);
        // renderTextureOne = new RenderTexture(1,1,24);
        renderTexture.enableRandomWrite = true;
        // renderTextureOne.enableRandomWrite = true;
        velocityGraphicsBuffer = new GraphicsBuffer(GraphicsBuffer.Target.Structured, nx*ny, sizeof(float)*2);
        vfx.SetGraphicsBuffer("VelocityBuffer",velocityGraphicsBuffer);
        vfx.SetInt("DIM_X",nx);
        vfx.SetInt("DIM_Y",ny);
        rhot = new ComputeBuffer((nx*ny+1)*2 + 2,sizeof(float));
        uv = new ComputeBuffer(nx*ny*2,sizeof(float));
        f = new ComputeBuffer(nx*ny*npop*2,sizeof(float));
        g = new ComputeBuffer(nx*ny*npop*2,sizeof(float));
        compute.SetInt("npop",npop);
        compute.SetInt("nx",nx);
        compute.SetInt("ny",ny);
        compute.SetInt("nxy9",nx*ny*9);
        compute.SetInt("nxy1",nx*ny + 1);
        // compute.SetFloat("aveRho",aveRho);
        compute.SetInt("radius",radius);
        compute.SetInt("initMode",(int)initMode);
        compute.SetFloat("initTr",initTr);
        OnValidate();

        init = compute.FindKernel("init");
        compute.SetBuffer(init,"rhot",rhot);
        compute.SetBuffer(init,"uv",uv);
        compute.SetBuffer(init,"f",f);
        compute.SetBuffer(init,"g",g);

        plotPass = compute.FindKernel("plotPass");
        compute.SetBuffer(plotPass,"g",g);
        // compute.SetBuffer(plotPass,"uv",uv);
        compute.SetBuffer(plotPass,"f",f);
        compute.SetBuffer(plotPass,"rhot",rhot);
        compute.SetTexture(plotPass,"renderTexture",renderTexture);

        plotSpeed = compute.FindKernel("plotSpeed");
        compute.SetBuffer(plotSpeed,"uv",uv);
        compute.SetTexture(plotSpeed,"renderTexture",renderTexture);

        collisionStreaming = compute.FindKernel("collisionStreaming");
        compute.SetBuffer(collisionStreaming,"rhot",rhot);
        compute.SetBuffer(collisionStreaming,"uv",uv);
        compute.SetBuffer(collisionStreaming,"velocityBuffer",velocityGraphicsBuffer);
        compute.SetBuffer(collisionStreaming,"f",f);

        collisionStreamingG = compute.FindKernel("collisionStreamingG");
        compute.SetBuffer(collisionStreamingG,"rhot",rhot);
        compute.SetBuffer(collisionStreamingG,"uv",uv);
        compute.SetBuffer(collisionStreamingG,"g",g);

        bouncebackBoundary = compute.FindKernel("bouncebackBoundary");
        compute.SetBuffer(bouncebackBoundary,"uv",uv);
        compute.SetBuffer(bouncebackBoundary,"rhot",rhot);
        compute.SetBuffer(bouncebackBoundary,"g",g);
        compute.SetBuffer(bouncebackBoundary,"f",f);

        oneLoop = compute.FindKernel("oneLoop");
        compute.SetBuffer(oneLoop,"rhot",rhot);

        compute.Dispatch(init,(nx*ny+63)/64,1,1);
    }

    private void FixedUpdate() 
    {
        for (int i = 0; i < loopCount; i++)
        {
            compute.Dispatch(plotPass,(nx*ny+63)/64,1,1);
            if(plotMode == PlotMode.Speed) compute.Dispatch(plotSpeed,(nx*ny+63)/64,1,1);
            compute.Dispatch(oneLoop,1,1,1);
            // if(grav!=0f)compute.Dispatch(oneLoop,1,1,1);
            compute.Dispatch(collisionStreaming,(nx+7)/8,(ny+7)/8,1);
            compute.Dispatch(collisionStreamingG,(nx+7)/8,(ny+7)/8,1);
            if(!periodic)compute.Dispatch(bouncebackBoundary,(nx+63)/64,1,1);
        }
        vfx.SetFloat("VelocityScale",((10f/(float)ny) * loopCount / Time.deltaTime) * Time.fixedDeltaTime/Time.deltaTime);
        RenderTexture.active = renderTexture;
        // renderTextureOne = renderTexture;
        // RenderTexture.active = renderTextureOne;
        plotTexture.ReadPixels(new Rect(0, 0, nx,ny), 0, 0);
        plotTexture.Apply();

        // Color[] pix = plotTexture.GetPixels();
        // print(pix[0].b == 0f);
    }
}
