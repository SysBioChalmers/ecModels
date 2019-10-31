pipeline {
  agent any
  stages {
    stage('Run') {
//      when {
//        branch 'develop'
//      }
      steps {
        withCredentials([string(credentialsId: 'GITHUB_TOKEN', variable: 'GITHUB_TOKEN')]) {
          sh '''
            python3 run.py
          '''
        }
      }
    }
  }
}
