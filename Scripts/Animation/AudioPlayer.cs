using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using System;

public class AudioPlayer : MonoBehaviour, IAnimationExecutor
{
    private AudioSource _audioSource;
    private float _startTime;
    private float? _endTime;

    public AudioPlayer(AudioSource audioSource, float startTime, float? endTime)
    {
        _audioSource = audioSource;
        _startTime = startTime;
        _endTime = endTime;
    }

    public void Execute()
    {
        _audioSource.Play();
    }

    public float GetStartTime()
    {
        return _startTime;
    }

    public float? GetEndTime()
    {
        return _endTime;
    }

}
