using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.EventSystems;
using UnityEngine.UI;
using System.Threading;

public class QuestionChoice : MonoBehaviour
{
    public List<GameObject> Questions;
    public List<GameObject> Answers;
    public List<GameObject> AnswersRight;

    private AudioSource audioGameObject;

    private GameObject activeQuestion;
    private GameObject activeAnswer;
    private GameObject activeQuestionOld;
    private GameObject activeAnswerOld;

    private int correctAnswers;

    // Use this for initialization
    void Start()
    {
        correctAnswers = 0;
        activeQuestionOld = null;
        activeQuestion = null;
        activeAnswerOld = null;
        activeAnswer = null;

        audioGameObject = transform.GetComponent<AudioSource>();
    }

    public void OnButtonClicked()
    {
        var activeObject = EventSystem.current.currentSelectedGameObject;

        if (activeObject.CompareTag("Question"))
        {
            ActivateGameObject(activeQuestion, activeObject);
            activeQuestion = activeObject;
        }
        else if (activeObject.CompareTag("Choice"))
        {
            ActivateGameObject(activeAnswer, activeObject);
            activeAnswer = activeObject;
        }

        if(activeQuestion != null && activeAnswer != null &&
            activeQuestion.GetComponent<Button>().enabled && activeAnswer.GetComponent<Button>().enabled)
        {
            var questionIndex = Questions.IndexOf(activeQuestion);
            var answerIndex = Answers.IndexOf(activeAnswer);
            var isCorrect = (AnswersRight[questionIndex] == activeAnswer);

            var line = transform.Find("Lines/" + IndexToName(questionIndex) + IndexToName(answerIndex)).gameObject;
            var lineAnimator = line.GetComponent<Animator>();

            lineAnimator.SetBool("IsLineIn", true);
            if (isCorrect)
            {
                correctAnswers++;
                lineAnimator.SetBool("IsLineCorrect", true);
                activeQuestion.GetComponent<Button>().enabled = false;
                activeAnswer.GetComponent<Button>().enabled = false;

                if(correctAnswers == 3)
                {
                    audioGameObject.PlayDelayed(6.0f);
                }
            }
            else
            {
                lineAnimator.SetBool("IsLineIncorrect", true);
            }

        }
    }

    private void ActivateGameObject(GameObject oldObject, GameObject activeObject)
    {
        if(oldObject != null)
        {
            var oldAnimator = oldObject.GetComponent<Animator>();
            oldAnimator.SetBool("IsOutlineIn", false);
            oldAnimator.SetBool("IsOutlineOut", true);
        }

        var activeAnimator = activeObject.GetComponent<Animator>();
        activeAnimator.SetBool("IsOutlineIn", true);
        activeAnimator.SetBool("IsOutlineOut", false);
    }

    string IndexToName(int index)
    {
        if(index == 0)
        {
            return "First";
        }else if(index == 1)
        {
            return "Second";
        }else if(index == 2)
        {
            return "Third";
        }
        return "";
    }
    
}
