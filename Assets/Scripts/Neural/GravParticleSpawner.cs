using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GravParticleSpawner : MonoBehaviour
{
    [SerializeField]
    GameObject gravParticlePrefab;
    [SerializeField]
    Transform startPoint;
    [SerializeField]
    Transform ansPoint;
    [SerializeField]
    Transform particleParent;
    [SerializeField]
    int particleCount;
    [SerializeField]
    float gravScale;
    [SerializeField]
    float dragScale;
    [SerializeField]
    float initVelScale;
    [SerializeField]
    int loopCount;
    [SerializeField]
    float thresDist;

    private void OnValidate() {
        GravParticle.GravScale = gravScale;
        GravParticle.LoopCount = loopCount;
        GravParticle.ThresDist = thresDist;
        GravParticle.DragScale = dragScale;
        GravParticle.InitVelScale = initVelScale;
    }

    void Start()
    {
        OnValidate();
        GravParticle.AnswerPoint = (RectTransform)ansPoint;
        for (int i = 0; i < particleCount; i++)
        {
            GameObject newParticle = Instantiate(gravParticlePrefab,particleParent);
            newParticle.transform.position = startPoint.position;
            GravParticle gp = newParticle.GetComponent<GravParticle>();
            float theta = Random.Range(0f,2f * Mathf.PI);
            gp.velocity = initVelScale * new Vector2(Mathf.Cos(theta),Mathf.Sin(theta));
            gp.initTheta = theta;
            SpriteRenderer sr = newParticle.GetComponent<SpriteRenderer>();
            sr.color = Random.ColorHSV(0f, 1f, 1f, 1f, 1f, 1f);
        }
    }

    void Update()
    {
        if(particleParent.childCount == 0)
        {
            TimeWriter.WriteTimes();
            UnityEditor.EditorApplication.isPlaying = false;
        }
    }
}
