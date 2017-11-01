using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.EventSystems;

[RequireComponent(typeof(AudioSource))]
public class AudioManager : MonoBehaviour, ITrackableSourceHandler, ITrackableStateHandler
{
    private bool isPlayed;
    private AudioSource audioSource;

    void Start()
    {
        isPlayed = false;
        audioSource = GetComponent<AudioSource>();
    }

    void Update()
    {

    }

    public void OnTrackableFound(GameObject gameObject)
    {
        if (audioSource)
        {
            if (audioSource.time > 0)
            {
                audioSource.UnPause();
            }
            else if(!isPlayed)
            {
                isPlayed = true;
                audioSource.Play();
            }
        }
    }

    public void OnTrackableLost(GameObject gameObject)
    {
        if (audioSource && audioSource.isPlaying)
            audioSource.Pause();
    }

    public void OnTrackableSourceChange()
    {
        /*if (audioSource)
        {
            isPlayed = false;
            audioSource.Stop();
        }*/
            
    }

}
