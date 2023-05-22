using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

public class DampedOccilator : MonoBehaviour
{
    public int LoopCount = 10;
    public float effZero;
    [SerializeField]
    Transform particle;
    [SerializeField]
    Transform origin;

    public float forceScale;
    public float dampScale;
    public float initPos;
    public float initVel;

    float vel;
    public float activeTime = 0f;
    bool done;

    public void Init()
    {
        done = false;
        activeTime = 0f;
        particle.position = new Vector3(0,initPos,0) + origin.position;
        vel = initVel;
    }

    void Start()
    {
        Init();
    }

    private void FixedUpdate() {
        if(done) return;
        for (int i = 0; i < LoopCount; i++)
        {
            float k1v = Force(particle.position.y, vel) * Time.fixedDeltaTime;
            float k1x = vel * Time.fixedDeltaTime;

            float k2v = Force(particle.position.y + k1x/2f,vel + k1v/2f) * Time.fixedDeltaTime;
            float k2x = (vel + k1v/2f) * Time.fixedDeltaTime;

            float k3v = Force(particle.position.y + k2x/2f,vel + k2v/2f) * Time.fixedDeltaTime;
            float k3x = (vel + k2v/2f) * Time.fixedDeltaTime;

            float k4v = Force(particle.position.y + k3x,vel + k3v) * Time.fixedDeltaTime;
            float k4x = (vel + k3v) * Time.fixedDeltaTime;

            vel += (k1v + 2f * k2v + 2f* k3v + k4v)/6f;
            particle.position += new Vector3(0,1,0) * (k1x + 2f * k2x + 2f* k3x + k4x)/6f;

            activeTime += Time.fixedDeltaTime;
            if(Done())
            {
                print(activeTime);
                done = true;
                break;
            }
        }
    }

    bool Done()
    {
        float pos = particle.position.y;
        float velocity = vel;
        if(pos*pos<=effZero*effZero&& velocity*velocity<=effZero*effZero)
        {
            return true;
        }
        return false;
    }

    float Force(float pos,float velocity)
    {
        return -forceScale*pos - dampScale *velocity; 
    }
}

[CustomEditor(typeof(DampedOccilator))]//拡張するクラスを指定
public class DampedOccilatorEditor : Editor {

  /// <summary>
  /// InspectorのGUIを更新
  /// </summary>
  public override void OnInspectorGUI(){
    //元のInspector部分を表示
    base.OnInspectorGUI ();
    DampedOccilator dampedOccilator = target as DampedOccilator;
    //ボタンを表示
    if (GUILayout.Button("Init")){
      dampedOccilator.Init();
    }  
  }

}  
