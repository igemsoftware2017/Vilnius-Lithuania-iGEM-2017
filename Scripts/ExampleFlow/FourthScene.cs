using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using System;
using UnityEngine.EventSystems;

public class FourthScene : MonoBehaviour, ITrackableStateHandler
{
    private List<IAnimationExecutor> animations;
    private float timeElapsed = 0.0f;
    private bool showAnimations = false;
    private int secondsCounter = 1;
    private float fadeSpeed = 0.1f;
    private List<string> madeAnswers;
    private int correctAnswers = 0;

    void Start()
    {
        madeAnswers = new List<string>();
        animations = new List<IAnimationExecutor>();
        var audioSource1 = transform.parent.Find("Audio").GetComponents<AudioSource>()[1];
        var audioSource2 = transform.parent.Find("Audio").GetComponents<AudioSource>()[2];
        var audioSource3 = transform.parent.Find("Audio").GetComponents<AudioSource>()[3];

        var playAudio1 = new AudioPlayer(audioSource1, 15.0f, null);
        var idle = new AnimatorParameterSetter<bool>(transform.GetComponent<Animator>(), "IsIdle", true, 14.0f, null);
        var shwoDna = new AnimatorParameterSetter<bool>(transform.GetComponent<Animator>(), "IsShowDna", true, 17.0f, null);
        var playAudio2 = new AudioPlayer(audioSource2, 26.0f, null);
        var correctDna = new AnimatorParameterSetter<bool>(transform.GetComponent<Animator>(), "IsCorrectDna", true, 26.0f, null);
        var bacteriaFill = new AnimatorParameterSetter<bool>(transform.GetComponent<Animator>(), "IsBacteria", true, 36.0f, null);
        var playAudio3 = new AudioPlayer(audioSource3, 41.0f, null);
        
        animations.Add(idle);
        animations.Add(playAudio1);
        animations.Add(shwoDna);
        animations.Add(playAudio2);
        animations.Add(correctDna);
        animations.Add(bacteriaFill);
        animations.Add(playAudio3);
    }

    void Update()
    {
        if (timeElapsed > 20.0 && timeElapsed < 25.0)
            return;

        if (!showAnimations)
            return;
        timeElapsed += Time.deltaTime;

        var animationsToUpdate = animations.Where(animation =>
            animation.GetStartTime() < timeElapsed &&
            (animation.GetEndTime() == null || animation.GetEndTime() > timeElapsed));

        if (animationsToUpdate.Count() == 0)
            return;

        animationsToUpdate.ToList().ForEach(animation => animation.Execute());

        var animationsToRemove = animations
            .Where(animation => animation.GetStartTime() < timeElapsed &&
            (animation.GetEndTime() == null || animation.GetEndTime() < timeElapsed));

        animationsToRemove.ToList().ForEach(animation => animations.Remove(animation));
    }


    public void OnTrackableFound(GameObject gameObject)
    {
        showAnimations = true;
    }

    public void OnTrackableLost(GameObject gameObject)
    {
        showAnimations = false;
    }

    public void OnButtonClicked(string ans)
    {
        if(madeAnswers.Contains(ans))
            return;

        if (ans.Contains("Correct"))
        {
            correctAnswers++;
        }

        madeAnswers.Add(ans);
        transform.Find("CorrectSigns/" + ans).gameObject.SetActive(true);

        if(correctAnswers == 2)
        {
            timeElapsed = 25.0f;
        }
    }
}
