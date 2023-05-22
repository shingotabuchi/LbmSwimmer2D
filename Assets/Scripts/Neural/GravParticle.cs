using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GravParticle : MonoBehaviour
{
    public static float GravScale;
    public static float LoopCount;
    public static float ThresDist;
    public static float DragScale;
    public static float InitVelScale;
    public static RectTransform AnswerPoint;
    public float drag;
    public Vector2 velocity;
    RectTransform rt;
    float spawnTime;
    public float initTheta;
    void Start()
    {
        spawnTime = Time.time;
        rt = GetComponent<RectTransform>();
    }

    void FixedUpdate()
    {
        for (int i = 0; i < LoopCount; i++)
        {
            Vector2 k1v = GravForce(rt.anchoredPosition,velocity) * Time.fixedDeltaTime;
            Vector2 k1x = velocity * Time.fixedDeltaTime;

            Vector2 k2v = GravForce(rt.anchoredPosition + k1x/2f,velocity + k1v/2f) * Time.fixedDeltaTime;
            Vector2 k2x = (velocity + k1v/2f) * Time.fixedDeltaTime;

            Vector2 k3v = GravForce(rt.anchoredPosition + k2x/2f,velocity + k2v/2f) * Time.fixedDeltaTime;
            Vector2 k3x = (velocity + k2v/2f) * Time.fixedDeltaTime;

            Vector2 k4v = GravForce(rt.anchoredPosition + k3x,velocity + k3v) * Time.fixedDeltaTime;
            Vector2 k4x = (velocity + k3v) * Time.fixedDeltaTime;

            velocity += (k1v + 2f * k2v + 2f* k3v + k4v)/6f;
            rt.anchoredPosition += (k1x + 2f * k2x + 2f* k3x + k4x)/6f;
            if(Done())
            {
                TimeWriter.times += (Time.time - spawnTime).ToString("0.000") + "\n";
                TimeWriter.thetas += (initTheta).ToString("0.000") + "\n";
                Destroy(gameObject);
                break;
            }
        }
    }

    Vector2 GravForce(Vector2 pos, Vector2 vel)
    {
        float dist = (AnswerPoint.anchoredPosition - pos).magnitude;
        return GravScale * (AnswerPoint.anchoredPosition - pos)/(dist*dist*dist) - DragScale*vel;
    }

    bool Done()
    {
        float sqrDist = (AnswerPoint.anchoredPosition - rt.anchoredPosition).sqrMagnitude;
        if(sqrDist < ThresDist*ThresDist) return true;
        return false;
    }
}
